##
#
# Single Rigid Body Traj Opt
#
##

# standard imports
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

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
        self.nu = 12   # [f_left, f_right, p_left, p_right]

    # SRB model continuous dynamics
    # https://arxiv.org/pdf/2207.04163
    def f_cont(self, x, u):

        # system parameters
        m = 35.0                                # mass [kg]
        g = 9.81                                # gravity [m/s^2]
        I = ca.diag(ca.DM([0.50, 0.50, 0.25]))  # body frame inertia matrix [kg*m^2]
        
        # extract the state
        p_com =  x[0:3]    # position in world frame
        quat =   x[3:7]    # orientation quaternion q_BW (body -> world), [w,x,y,z]
        v_com =  x[7:10]   # linear velocity in world frame
        w_body = x[10:13]  # body frame angular velocity

        # extract the inputs (world frame)
        f_left =  u[0:3]  # left foot force
        f_right = u[3:6]  # right foot force
        p_left =  u[6:9]  # left foot position
        p_right = u[9:12] # right foot position

        # rotation of body expressed in world frame
        R_BW = self._quat_to_rotmat(quat)

        # net force in the world frame
        F_net_W = f_left + f_right + ca.DM([0, 0, -m*g])

        # net moment about COM
        M_net_W = ca.cross(p_left - p_com, f_left) + ca.cross(p_right - p_com, f_right)
        M_net_B = R_BW.T @ M_net_W  # express moment in body frame

        # translational dynamics
        p_com_dot = v_com
        v_com_dot = (1.0 / m) * F_net_W

        # quaternion rate
        w_body_quat = ca.vertcat(0, w_body)  # augment angular velocity to quaternion form [0, wx, wy, wz]
        quat_dot = 0.5 * self._quat_hamilton(quat, w_body_quat)
        
        # angular dynamics
        w_body_dot = ca.solve(I, M_net_B - ca.cross(w_body, I @ w_body))

        # build the dynamics vector
        x_dot = ca.vertcat(
            p_com_dot,
            quat_dot,
            v_com_dot,
            w_body_dot
        )

        return x_dot
    
    # SRB model discrete dynamics using Euler integration
    def f_disc(self, x, u, dt):
        
        # get the continuous dynamics vector
        x_dot = self.f_cont(x, u)

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

    # running cost
    def running_cost(self, x, u, x_goal):

        # simple quadratic 
        w_pos = 10.0
        w_ori = 10.0
        w_vel = 1.0
        w_omega = 1.0
        Qx = ca.diag(ca.vertcat(
            w_pos, w_pos, w_pos,        # p_com
            w_ori, w_ori, w_ori,        # quat vector part
            w_vel, w_vel, w_vel,        # v_com
            w_omega, w_omega, w_omega   # w_body
        ))

        # compute errors
        pos_err = x[0:3] - x_goal[0:3]
        vel_err = x[7:10] - x_goal[7:10]
        omega_err = x[10:13] - x_goal[10:13]
        quat_err = self._quat_error(x[3:7], x_goal[3:7])
        quat_err = ca.if_else(quat_err[0] < 0, -quat_err, quat_err)

        # state error vector
        e_x = ca.vertcat(
            pos_err,
            quat_err[1:4],  # vector part of quaternion error
            vel_err,
            omega_err
        )

        # state cost
        cost_state = e_x.T @ Qx @ e_x
        
        # penalize forces and positions
        w_force = 0.1
        w_feet = 0.0
        Qu = ca.diag(ca.vertcat(
            w_force, w_force, w_force,    # f_left
            w_force, w_force, w_force,    # f_right
            w_feet, w_feet, w_feet,       # p_left
            w_feet, w_feet, w_feet        # p_right
        ))

        # compute errors
        f_left = u[0:3]
        f_right = u[3:6]
        p_left_W = u[6:9]
        p_right_W = u[9:12]

        # p_left_B = p_left - x[0:3]
        # p_right_B = p_right - x[0:3]

        # input error vector
        e_u = ca.vertcat(
            f_left,
            f_right,
            p_left_W,
            p_right_W
        )

        # input cost
        cost_input = e_u.T @ Qu @ e_u

        # total cost
        cost_tot = cost_input + cost_state

        return cost_tot
    
    # terminal cost
    def terminal_cost(self, x, x_goal):

        # terminal weights (usually larger than running)
        w_pos   = 200.0
        w_ori   = 200.0
        w_vel   = 50.0
        w_omega = 50.0

        QT = ca.diag(ca.vertcat(
            w_pos, w_pos, w_pos,        # p_com
            w_ori, w_ori, w_ori,        # quat vector part
            w_vel, w_vel, w_vel,        # v_com
            w_omega, w_omega, w_omega   # w_body
        ))

        # compute errors (same as running_cost)
        pos_err   = x[0:3]   - x_goal[0:3]
        vel_err   = x[7:10]  - x_goal[7:10]
        omega_err = x[10:13] - x_goal[10:13]

        quat_err = self._quat_error(x[3:7], x_goal[3:7])
        quat_err = ca.if_else(quat_err[0] < 0, -quat_err, quat_err)

        eT = ca.vertcat(
            pos_err,
            quat_err[1:4],   # vector part of quaternion error
            vel_err,
            omega_err
        )

        cost_terminal = eT.T @ QT @ eT
        return cost_terminal

##############################################################

# index helpers (readability)
IDX_P      = slice(0, 3)
IDX_Q      = slice(3, 7)
IDX_V      = slice(7, 10)
IDX_W      = slice(10, 13)

IDX_FL_X   = 0
IDX_FL_Y   = 1
IDX_FL_Z   = 2
IDX_FR_X   = 3
IDX_FR_Y   = 4
IDX_FR_Z   = 5
IDX_FL     = slice(0, 3)
IDX_FR     = slice(3, 6)

IDX_PL     = slice(6, 9)
IDX_PR     = slice(9, 12)

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
dt = 0.04        # time step
T = 4.0          # total time
N = int(T / dt)  # number of intervals

# make the optimizer
opti = ca.Opti()

# horizon variables
X = opti.variable(nx, N + 1)  # states over the horizon
U = opti.variable(nu, N)      # inputs over the horizon

# initial condition
x0 = np.array([0, 0, 1.0,  # p_com
               1, 0, 0, 0, # quaternion
               0, 0, 0,    # v_com
               0, 0, 0])   # w_body

# desired goal state
pitch_goal = np.deg2rad(30) 
x_goal = np.array([1.0, 0, 0.5, # p_com
                   np.cos(pitch_goal/2), 0, np.sin(pitch_goal/2), 0,  # quaternion
                   0, 0, 0,     # v_com
                   0, 0, 0])    # w_body

# set the initial condition 
opti.subject_to(X[:, 0] == x0)

# system dynamics constraints at each time step
for k in range(N):
    x_next = f(X[:, k], U[:, k], dt)
    opti.subject_to(X[:, k + 1] == x_next)

# input bounds
f_max = 1000.0    
opti.subject_to(opti.bounded(-f_max, U[IDX_FL, :], f_max))    # f_left
opti.subject_to(opti.bounded(-f_max, U[IDX_FR, :], f_max))    # f_right
opti.subject_to(U[IDX_FL_Z, :] >= 0)  # f_left_z
opti.subject_to(U[IDX_FR_Z, :] >= 0)  # f_right_z

mu = 1.0
fxL, fyL, fzL = U[IDX_FL_X, :], U[IDX_FL_Y, :], U[IDX_FL_Z, :]
fxR, fyR, fzR = U[IDX_FR_X, :], U[IDX_FR_Y, :], U[IDX_FR_Z, :]
opti.subject_to(fxL**2 + fyL**2 <= (mu*fzL)**2)
opti.subject_to(fxR**2 + fyR**2 <= (mu*fzR)**2)

# foot positions: fix them (then you don't need bounds/cost on them)
p_left0  = ca.DM([0.1,  0.1, 0.0])
p_right0 = ca.DM([0.1, -0.1, 0.0])
opti.subject_to(U[IDX_PL, :] == ca.repmat(p_left0,  1, N))
opti.subject_to(U[IDX_PR, :] == ca.repmat(p_right0, 1, N))

# objective function 
J = 0
for k in range(N):
    J += srb.running_cost(X[:, k], U[:, k], x_goal)
J += srb.terminal_cost(X[:, N], x_goal)
opti.minimize(J)

# -----------------------------
# Initial guesses (highly recommended)
# -----------------------------
opti.set_initial(X, np.tile(x0.reshape(-1, 1), (1, N+1)))
opti.set_initial(U, 0)

# better force guess: support weight evenly
m = 35.0
g = 9.81
opti.set_initial(U[2, :], 0.5*m*g)  # fLz
opti.set_initial(U[5, :], 0.5*m*g)  # fRz

opti.solver("ipopt", {"expand": True}, {"max_iter": 3000, "tol": 1e-6})
sol = opti.solve()

X_sol = sol.value(X)
U_sol = sol.value(U)

# plot some results
plt.figure()
plt.subplot(3,1,1)
plt.plot(X_sol[0, :], label='x')
plt.plot(X_sol[1, :], label='y')
plt.plot(X_sol[2, :], label='z')
plt.title('CoM Position')
plt.legend()

plt.subplot(3,1,2)
plt.plot(X_sol[7, :], label='vx')
plt.plot(X_sol[8, :], label='vy')
plt.plot(X_sol[9, :], label='vz')
plt.title('CoM Velocity')
plt.legend()

plt.show()