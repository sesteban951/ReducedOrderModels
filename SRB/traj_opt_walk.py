##
#
# Single Rigid Body Traj Opt
#
##

# standard imports
import numpy as np
import matplotlib.pyplot as plt
import os

# casadi import
import casadi as ca

##############################################################
# Single Rigid Body Dynamics
##############################################################

class SRBDynamics:

    # initialize the class
    def __init__(self):
        
        # state dimension
        self.nq = 7  # [p_com, quat]
        self.nv = 6  # [v_com, w_body]

        # input dimension
        self.nu = 12   # [F_left, F_right, M_left, M_right]

        # system parameters
        self.m = 35.0                                # G1 mass [kg]
        self.g = 9.81                                # gravity [m/s^2]
        self.I = ca.diag(ca.DM([0.50, 0.50, 0.25]))  # body frame inertia matrix [kg*m^2]
        self.hip_offset = 0.1185                     # G1 hip offset from Base [m]

        # simple quadratic 
        w_pos = 10.0
        w_ori = 8.0
        w_vel = 2.0
        w_omega = 2.0
        self.Qx = ca.diag(ca.vertcat(
            w_pos, w_pos, w_pos,        # p_com
            w_ori, w_ori, w_ori,        # quat vector part
            w_vel, w_vel, w_vel,        # v_com
            w_omega, w_omega, w_omega   # w_body
        ))

        # penalize forces and positions
        w_force = 0.01
        w_moment = 0.01
        self.Qu = ca.diag(ca.vertcat(
            w_force, w_force, w_force,    # F_left
            w_force, w_force, w_force,    # F_right
            w_moment, w_moment, w_moment, # M_left
            w_moment, w_moment, w_moment  # M_right
        ))

        # terminal weights (usually larger than running)
        self.Qx_f = 500.0 * self.Qx

    ###############################################################
    # Dynamics
    ###############################################################

    # SRB model continuous dynamics
    # https://arxiv.org/pdf/2207.04163
    def f_cont(self, x, u, p_left_W, p_right_W):
        
        # extract the state
        p_com =  x[0:3]    # position in world frame
        quat =   x[3:7]    # orientation quaternion q_BW (body -> world), [w,x,y,z]
        v_com =  x[7:10]   # linear velocity in world frame
        w_body = x[10:13]  # body frame angular velocity

        # extract the inputs (world frame)
        F_left =  u[0:3]  # left foot force    in world frame
        F_right = u[3:6]  # right foot force   in world frame
        M_left =  u[6:9]  # left foot moment   in world frame
        M_right = u[9:12] # right foot moment  in world frame

        # rotation of body expressed in world frame
        R_BW = self._quat_to_rotmat(quat)

        # net force in the world frame
        F_net_W = F_left + F_right + ca.DM([0, 0, -self.m * self.g])

        # net moment about COM
        M_net_W = (ca.cross(p_left_W - p_com, F_left) 
                 + ca.cross(p_right_W - p_com, F_right)
                 + M_left + M_right)
        M_net_B = R_BW.T @ M_net_W  # express moment in body frame

        # translational dynamics
        p_com_dot = v_com
        v_com_dot = (1.0 / self.m) * F_net_W

        # quaternion rate
        w_body_quat = ca.vertcat(0, w_body)  # augment angular velocity to quaternion form [0, wx, wy, wz]
        quat_dot = 0.5 * self._quat_hamilton(quat, w_body_quat)
        
        # angular dynamics
        w_body_dot = ca.solve(self.I, M_net_B - ca.cross(w_body, self.I @ w_body))

        # build the dynamics vector
        x_dot = ca.vertcat(
            p_com_dot,
            quat_dot,
            v_com_dot,
            w_body_dot
        )

        return x_dot
    
    # SRB model discrete dynamics using Euler integration
    def f_disc(self, x, u, dt, p_left_W, p_right_W):
        
        # get the continuous dynamics vector
        x_dot = self.f_cont(x, u, p_left_W, p_right_W)

        # Euler integration
        x_next = x + dt * x_dot

        # normalize the quaternion
        quat_next = x_next[3:7]
        quat_next = quat_next / (ca.norm_2(quat_next) + 1e-12)
        
        # rebuild state to avoid slice assignment
        x_next = ca.vertcat(
            x_next[0:3],     # p
            quat_next,       # quat
            x_next[7:10],    # v
            x_next[10:13]    # w_body
        )
        
        return x_next
    
    ###############################################################
    # Helper Functions
    ###############################################################

    # quaternion to rotation matrix
    def _quat_to_rotmat(self, q):

        # extract quaternion components
        qw = q[0]
        qx = q[1]
        qy = q[2]
        qz = q[3]
        
        # compute rotation matrix
        R_00 = 1 - 2*(qy*qy + qz*qz)
        R_01 = 2*(qx*qy - qz*qw)
        R_02 = 2*(qx*qz + qy*qw)
        R_10 = 2*(qx*qy + qz*qw)
        R_11 = 1 - 2*(qx*qx + qz*qz)
        R_12 = 2*(qy*qz - qx*qw)
        R_20 = 2*(qx*qz - qy*qw)
        R_21 = 2*(qy*qz + qx*qw)
        R_22 = 1 - 2*(qx*qx + qy*qy)
        R = ca.vertcat(
            ca.horzcat(R_00, R_01, R_02),
            ca.horzcat(R_10, R_11, R_12),
            ca.horzcat(R_20, R_21, R_22)
        )

        return R
    
    # rotation matrix to quaternion
    def _rotmat_to_euler_ZYX(self, R):

        R_20 = R[2, 0]
        R_20 = np.clip(R_20, -1.0, 1.0)  # numerical safety

        # compute pitch
        p = np.arcsin(-R_20)

        # check gimbal lock
        eps = 1e-8
        if abs(np.cos(p)) > eps:
            r = np.arctan2(R[2, 1], R[2, 2])
            y  = np.arctan2(R[1, 0], R[0, 0])
        else:
            # gimbal lock: yaw-roll coupling
            # set roll = 0 and compute yaw from other terms
            r = 0.0
            y  = np.arctan2(-R[0, 1], R[1, 1])

        return y, p, r

    
    # quaternion Hamilton product
    def _quat_hamilton(self, a, b):
        """
        Hamilton product c = a ⊗ b for quaternions in [w, x, y, z] convention.
        Returns a 4x1 SX vector.
        """

        # unpck the quaternions
        aw, ax, ay, az = a[0], a[1], a[2], a[3]
        bw, bx, by, bz = b[0], b[1], b[2], b[3]

        # compute the Hamilton product
        c = ca.vertcat(
            aw*bw - ax*bx - ay*by - az*bz,
            aw*bx + ax*bw + ay*bz - az*by,
            aw*by - ax*bz + ay*bw + az*bx,
            aw*bz + ax*by - ay*bx + az*bw
        )

        return c
    
    # quaternion conjugate
    def _quat_conj(self, q):
        return ca.vertcat(q[0], -q[1], -q[2], -q[3])

    # quaternion error
    def _quat_error(self, a, b):
        
        # compute the conjugate of a
        a_conj = self._quat_conj(a)

        # compute the error quaternion
        q_err = self._quat_hamilton(b, a_conj)

        return q_err
    
    # friction cone matrix for single force
    def _friction_cone_matrix(self, mu):

        # build the friction cone matrix
        A = ca.vertcat(
            ca.horzcat( 1,  0, -mu),
            ca.horzcat(-1,  0, -mu),
            ca.horzcat( 0,  1, -mu),
            ca.horzcat( 0, -1, -mu),
            ca.horzcat( 0,  0,  -1)
        )
        b = ca.DM.zeros(5, 1)

        return A, b
    
    # build a contact schedule
    def contact_schedule(self, dt, N, T_SSP, init_stance="L"):

        # compute the number of nodes per single support phase
        N_node_SSP = int(round(T_SSP / dt))
        N_node_SSP = max(N_node_SSP, 1)  # at least 1 step

        # storage for the contact schedule
        contact_L = np.zeros(N, dtype=int)
        contact_R = np.zeros(N, dtype=int)
        phase_idx = np.zeros(N, dtype=int)

        # initialize the stance
        stance = init_stance.upper()
        assert stance in ["L", "R"], "init_stance must be 'L' or 'R'"

        # loop through trajectory length
        for k in range(N):

            # phase number
            phase = k // N_node_SSP
            phase_idx[k] = phase

            if stance == "L":
                contact_L[k] = 1
            else:
                contact_R[k] = 1

            # flip at the start of each new phase
            if (k+1) % N_node_SSP == 0:
                stance = "R" if stance == "L" else "L"

        return contact_L, contact_R, phase_idx, N_node_SSP

    # get the actual contact state
    def get_contact_state(self, t, T_SSP, init_stance="L"):

        # clamp negative time
        if t < 0:
            t = 0.0

        # compute the phase
        phase = int(np.floor(t / T_SSP))

        # determine which foot is in stance
        left_stance = (init_stance.upper() == "L")
        if phase % 2 == 1:
            left_stance = not left_stance

        # build contact state
        cL = 1 if left_stance else 0
        cR = 0 if left_stance else 1

        return cL, cR, phase

    # get foot positions at phase k
    def p_left_W_at(self, k, phase_idx):
        # get which current step we are in
        curr_step = int(phase_idx[k])

        # return the foot position
        p_left_W_k = ca.vertcat(P_L_xy[:, curr_step], 0.0)

        return p_left_W_k

    def p_right_W_at(self, k, phase_idx):
        # get which current step we are in
        curr_step = int(phase_idx[k])

        # return the foot position
        p_right_W_k = ca.vertcat(P_R_xy[:, curr_step], 0.0)

        return p_right_W_k

    ###############################################################
    # Cost Functions
    ###############################################################

    # running cost
    def running_cost(self, x, u, x_goal):

        # compute errors
        pos_err = x[0:3] - x_goal[0:3]
        vel_err = x[7:10] - x_goal[7:10]
        omega_err = x[10:13] - x_goal[10:13]

        quat_err = self._quat_error(x[3:7], x_goal[3:7])
        quat_err = ca.if_else(quat_err[0] < 0, -quat_err, quat_err)
        quat_err_img = quat_err[1:4]

        # state error vector
        e_x = ca.vertcat(
            pos_err,
            quat_err_img,
            vel_err,
            omega_err
        )

        # state cost
        cost_state = 0.5 * e_x.T @ self.Qx @ e_x

        # compute errors
        F_left = u[0:3]
        F_right = u[3:6]
        M_left_W = u[6:9]
        M_right_W = u[9:12]

        # input error vector
        e_u = ca.vertcat(
            F_left,
            F_right,
            M_left_W,
            M_right_W
        )

        # input cost
        cost_input = 0.5 * e_u.T @ self.Qu @ e_u

        # total cost
        cost_tot = cost_input + cost_state

        return cost_tot
    
    # terminal cost
    def terminal_cost(self, x, x_goal):

        # compute errors (same as running_cost)
        pos_err   = x[0:3]   - x_goal[0:3]
        vel_err   = x[7:10]  - x_goal[7:10]
        omega_err = x[10:13] - x_goal[10:13]

        quat_err = self._quat_error(x[3:7], x_goal[3:7])
        quat_err = ca.if_else(quat_err[0] < 0, -quat_err, quat_err)
        quat_err_img = quat_err[1:4]
        
        # state error vector
        e = ca.vertcat(
            pos_err,
            quat_err_img,
            vel_err,
            omega_err
        )

        # terminal cost
        cost_terminal = e.T @ self.Qx_f @ e

        return cost_terminal

##############################################################

# state input indices
IDX_P      = slice(0, 3)
IDX_Q      = slice(3, 7)
IDX_V      = slice(7, 10)
IDX_W      = slice(10, 13)

# force input indices
IDX_FL_X   = 0
IDX_FL_Y   = 1
IDX_FL_Z   = 2
IDX_FR_X   = 3
IDX_FR_Y   = 4
IDX_FR_Z   = 5
IDX_FL     = slice(0, 3)
IDX_FR     = slice(3, 6)

# moment input indices
IDX_ML_X   = 6
IDX_ML_Y   = 7
IDX_ML_Z   = 8
IDX_MR_X   = 9
IDX_MR_Y   = 10
IDX_MR_Z   = 11
IDX_ML     = slice(6, 9)
IDX_MR     = slice(9, 12)

##############################################################
# Trajectory Optimization
##############################################################

# create the dynamics object
srb = SRBDynamics()
f = srb.f_disc
nq = srb.nq
nv = srb.nv
nx = nq + nv
nu = srb.nu

# optimization settings
dt = 0.025       # time step
T = 5.0          # total time
N = int(T / dt)  # number of intervals

# create the contact schedule
T_SSP = 0.4
schedule = srb.contact_schedule(dt, N, T_SSP, init_stance="L")
contact_L, contact_R, phase_idx, N_node_SSP = schedule
N_steps = int(phase_idx[-1] + 1)  # number of SSP phases actually used

# ----------------------------------------------------------
# Setup the optimization problem
# ----------------------------------------------------------

# make the optimizer
opti = ca.Opti()

# horizon variables
X = opti.variable(nx, N + 1)        # states over the horizon
U = opti.variable(nu, N)            # inputs over the horizon
P_L_xy = opti.variable(2, N_steps)  # left foot (x,y) for each SSP phase
P_R_xy = opti.variable(2, N_steps)  # right foot (x,y) for each SSP phase

# initial condition
x0 = np.array([0, 0, 1,    # p_com
               1, 0, 0, 0, # quaternion
               0, 0, 0,    # v_com
               0, 0, 0])   # w_body

# desired goal state
x_goal = np.array([1.0, 1.0, 1,  # p_com
                   1, 0, 0, 0,  # quaternion
                   0, 0, 0,     # v_com
                   0, 0, 0])    # w_body

# set the initial condition 
opti.subject_to(X[:, 0] == x0)

# system dynamics constraints at each time step
for k in range(N):
    x_next = srb.f_disc(X[:, k], U[:, k], dt, srb.p_left_W_at(k, phase_idx), srb.p_right_W_at(k, phase_idx))
    opti.subject_to(X[:, k + 1] == x_next)

# foot position constraints across phase boundaries
for s in range(1, N_steps):

    k0 = s * N_node_SSP

    if k0 >= N:
        break

    if contact_L[k0] == 1:
        # left is stance in phase s -> left cannot change across boundary
        opti.subject_to(P_L_xy[:, s] == P_L_xy[:, s-1])
        # right is swing -> allowed to change
    else:
        # right is stance in phase s
        opti.subject_to(P_R_xy[:, s] == P_R_xy[:, s-1])
        # left is swing -> allowed to change

# state constraints
z_min = 0.2
for k in range(N+1):
    opti.subject_to(X[2, k] >= z_min)  # z com min height

# force limits
m_max = 500.0  
mu = 1.0
A, b = srb._friction_cone_matrix(mu)
for k in range(N):
    if contact_L[k] == 1:
        opti.subject_to(A @ U[IDX_FL, k] <= b)
        opti.subject_to(opti.bounded(-m_max, U[IDX_ML, k], m_max))
    else:
        opti.subject_to(U[IDX_FL, k] == 0)
        opti.subject_to(U[IDX_ML, k] == 0)

    if contact_R[k] == 1:
        opti.subject_to(A @ U[IDX_FR, k] <= b)
        opti.subject_to(opti.bounded(-m_max, U[IDX_MR, k], m_max))
    else:
        opti.subject_to(U[IDX_FR, k] == 0)
        opti.subject_to(U[IDX_MR, k] == 0)

# objective function 
J = 0
for k in range(N):
    J += srb.running_cost(X[:, k], U[:, k], x_goal)

# either terminal cost or final state constraint
# J += srb.terminal_cost(X[:, N], x_goal)
opti.subject_to(X[:, N] == x_goal)

# set the objective
opti.minimize(J)

opti.set_initial(X, np.tile(x0.reshape(-1, 1), (1, N+1)))
opti.set_initial(U, 0)

# better force guess: support weight evenly
opti.set_initial(U[2, :], 0.5 * srb.m * srb.g)  # fLz
opti.set_initial(U[5, :], 0.5 * srb.m * srb.g)  # fRz

opti.set_initial(P_L_xy, np.tile(np.array([[0.0],[ srb.hip_offset]]), (1, N_steps)))
opti.set_initial(P_R_xy, np.tile(np.array([[0.0],[-srb.hip_offset]]), (1, N_steps)))

# ----------------------------------------------------------
# Solve the optimization
# ----------------------------------------------------------

# solver settings
opti.solver(
    "ipopt",
)
sol = opti.solve()
X_sol = sol.value(X)
U_sol = sol.value(U)
P_L_xy_sol = sol.value(P_L_xy)
P_R_xy_sol = sol.value(P_R_xy)

# ----------------------------------------------------------
# Save
# ----------------------------------------------------------

# create the time array
time = np.linspace(0, T, N+1)

BIG = 1e6

P_L_W = BIG * np.ones((3, N + 1))
P_R_W = BIG * np.ones((3, N + 1))
P_support_W = BIG * np.ones((3, N + 1))

for k in range(N + 1):
    # phase index: defined on intervals (0..N-1)
    s = int(phase_idx[k]) if k < len(phase_idx) else int(phase_idx[-1])

    # contacts: also defined on intervals (0..N-1)
    cL = int(contact_L[k]) if k < len(contact_L) else int(contact_L[-1])
    cR = int(contact_R[k]) if k < len(contact_R) else int(contact_R[-1])

    pL_xy = P_L_xy_sol[:, s]
    pR_xy = P_R_xy_sol[:, s]

    if cL == 1:
        P_L_W[:, k] = np.array([pL_xy[0], pL_xy[1], 0.0])
        P_support_W[:, k] = P_L_W[:, k]
    else:
        P_L_W[:, k] = np.array([BIG, BIG, BIG])

    if cR == 1:
        P_R_W[:, k] = np.array([pR_xy[0], pR_xy[1], 0.0])
        if cL == 0:
            P_support_W[:, k] = P_R_W[:, k]
    else:
        P_R_W[:, k] = np.array([BIG, BIG, BIG])

# (N+1,3) for saving/plotting
P_L_W_T = P_L_W.T
P_R_W_T = P_R_W.T
P_sup_T = P_support_W.T


# save the solution as csv
X_sol_T = X_sol.T
U_sol_T = U_sol.T
P_L_xy_sol_T = P_L_xy_sol.T
P_R_xy_sol_T = P_R_xy_sol.T
save_dir = "./SRB/results/walk/"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
time_file =  "./SRB/results/walk/time.csv"
state_file = "./SRB/results/walk/states.csv"
input_file = "./SRB/results/walk/inputs.csv"
p_left_file = "./SRB/results/walk/p_left.csv"
p_right_file = "./SRB/results/walk/p_right.csv"
p_sup_file = "./SRB/results/walk/p_support.csv"
np.savetxt(time_file, time, delimiter=",")
np.savetxt(state_file, X_sol_T, delimiter=",")
np.savetxt(input_file, U_sol_T, delimiter=",")
np.savetxt(p_left_file, P_L_W_T, delimiter=",")
np.savetxt(p_right_file, P_R_W_T, delimiter=",")
np.savetxt(p_sup_file, P_sup_T, delimiter=",")
print(f"Saved time to {time_file}")
print(f"Saved states to {state_file}")
print(f"Saved inputs to {input_file}")
print(f"Saved left foot positions to {p_left_file}")
print(f"Saved right foot positions to {p_right_file}")
print(f"Saved support foot positions to {p_sup_file}")