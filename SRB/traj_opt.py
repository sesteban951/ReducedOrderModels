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
    p_com = x[0:3]    # position
    theta = x[3:6]    # euler angles
    v_com = x[6:9]    # linear velocity
    w_body = x[9:12]  # angular velocity

    # extract the inputs
    f_left =  u[0:3]  # control torques
    f_right = u[3:6]  # control torques

    # extract the euler angles (world frame)
    r = theta[0]
    p = theta[1]
    y = theta[2]

    # create the rotation matrix from body to world frame
    R = euler_to_rotmat(r, p, y)

    # precompute sines and cosines
    cr, sr = ca.cos(r), ca.sin(r)
    cp, sp = ca.cos(p), ca.sin(p)














##############################################################


