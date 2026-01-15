##
#
# Single Rigid Body Traj Opt
#
##

# standard imports
import numpy as np
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
        F_net_W = f_left + f_right + ca.SX([0, 0, -m*g])

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

    # quaternion to rotation matrix
    def _quat_to_rotmat(self, q):

        # extract quaternion components
        qw = q[0]
        qx = q[1]
        qy = q[2]
        qz = q[3]
        
        # compute rotation matrix
        R = ca.SX.zeros(3, 3)
        R[0, 0] = 1 - 2*(qy**2 + qz**2)
        R[0, 1] = 2*(qx*qy - qz*qw)
        R[0, 2] = 2*(qx*qz + qy*qw)
        R[1, 0] = 2*(qx*qy + qz*qw)
        R[1, 1] = 1 - 2*(qx**2 + qz**2)
        R[1, 2] = 2*(qy*qz - qx*qw)
        R[2, 0] = 2*(qx*qz - qy*qw)
        R[2, 1] = 2*(qy*qz + qx*qw)
        R[2, 2] = 1 - 2*(qx**2 + qy**2)

        return R

##############################################################
# Cost Functions
##############################################################

# running cost
def running_cost(x, u, x_goal):

    # simple quadratic 
    w_pos = 10.0
    w_ori = 10.0
    w_vel = 1.0
    w_omega = 1.0
    Qx = ca.diag(ca.vertcat(
        w_pos, w_pos, w_pos,        # p_com
        w_ori, w_ori, w_ori, w_ori, # quaternion # TODO: consider orientation error differently
        w_vel, w_vel, w_vel,        # v_com
        w_omega, w_omega, w_omega   # w_body
    ))
    
    # penalize forces and positions
    w_force = 0.1
    







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
dt = 0.01        # time step
T = 5.0          # total time
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
x_goal = np.array([1.0, 0, 0.5, # p_com
                   0.707, 0, 0.707, 0,  # quaternion
                   0, 0, 0,     # v_com
                   0, 0, 0])    # w_body

# set the initial condition 
opti.subject_to(X[:, 0] == x0)

# set the final condition
opti.subject_to(X[:, N] == x_goal)

# system dynamics constraints at each time step
for k in range(N):
    x_next = f(X[:, k], U[:, k], dt)
    opti.subject_to(X[:, k + 1] == x_next)


