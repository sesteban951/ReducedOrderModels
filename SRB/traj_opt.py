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

# skew symmetric function
def skew(v):

    # return skew symmetric matrix
    v1 = v[0]; v2 = v[1]; v3 = v[2]
    S = ca.vertcat(
        ca.horzcat(  0, -v3,  v2),
        ca.horzcat( v3,   0, -v1),
        ca.horzcat(-v2,  v1,   0),
    )

    return S

# compute the rotation matrix from euler angles
def euler_to_rotmat(r, p, y):

    # precompute the sines and cosines
    cy, sy = ca.cos(y), ca.sin(y)
    cp, sp = ca.cos(p), ca.sin(p)
    cr, sr = ca.cos(r), ca.sin(r)

    # compute the rotation matrix
    Rz = ca.vertcat(
        ca.horzcat( cy, -sy, 0),
        ca.horzcat( sy,  cy, 0),
        ca.horzcat(  0,   0, 1),
    )
    Ry = ca.vertcat(
        ca.horzcat( cp, 0, sp),
        ca.horzcat(  0, 1,  0),
        ca.horzcat(-sp, 0, cp),
    )
    Rx = ca.vertcat(
        ca.horzcat( 1,   0,   0),
        ca.horzcat( 0,  cr, -sr),
        ca.horzcat( 0,  sr,  cr),
    )

    # compute the overall rotation matrix
    R = Rz @ Ry @ Rx

    return R

# dynamics function
def f_cont(self, x, u):

    # system parameters
    m = 35.0   # mass [kg]
    g = 9.81   # gravity [m/s^2]
    I = ca.diag(ca.DM([0.25, 0.25, 0.50]))  # inertia matrix [kg*m^2]
    
    # extract the state
    p_com = x[0:3]    # position in world frame
    theta = x[3:6]    # euler angles (roll, pitch, yaw) in world frame
    v_com = x[6:9]    # linear velocity in world frame
    w_body = x[9:12]  # angular velocity in body frame

    # extract the inputs
    f_left =  u[0:3]  # control forces  in world frame
    f_right = u[3:6]  # control torques in world frame

    # extract the euler angles (world frame)
    r = theta[0]  # roll
    p = theta[1]  # pitch
    y = theta[2]  # yaw

    # create the rotation matrix from body to world frame
    R = euler_to_rotmat(r, p, y) 

    # body angular velocity to euler angular rates
    cr, sr = ca.cos(r), ca.sin(r)
    cp, sp = ca.cos(p), ca.sin(p)
    T = ca.vertcat(
        ca.horzcat(1, sr*sp/cp, cr*sp/cp),
        ca.horzcat(0,       cr,      -sr),
        ca.horzcat(0,    sr/cp,    cr/cp),
    )
    theta_dot = T @ w_body

    # net force
    p_left = ca.DM([-0.5, 0.0, 0.0])   # left foot position in world frame
    p_right = ca.DM([0.5, 0.0, 0.0])   # right foot position in world frame

    F_net = f_left + f_right
    M_net = skew(p_left - p_com) @ f_left + skew(p_right - p_com) @ f_right

    # translational dynamics
    g_vec = ca.vertcat([0.0, 0.0, -g])
    p_dot = v_com
    v_dot = (1/m) * F_net + g_vec

    # rotational dynamics
    M_body = R.T @ M_net

    # w_dot = I^{-1} * (M_body - w_body x (I * w_body))
    w_dot = ca.solve(I, (M_body - skew(w_body)) @ (I @ w_body))

    # build dynamics vector
    x_dot = ca.vertcat(
        p_dot,
        theta_dot,
        v_dot,
        w_dot,
    )

    return x_dot











##############################################################


