%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC implementation of HLIP using YALMIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
yalmip('clear');

% HLIP Parameters
T_DSP = 0;
T_SSP = 0.15;
g = 9.81;
z0 = 0.5;
lambda = sqrt(g/z0);
A_SSP = [0, 1;
         lambda^2, 0];

% step-to-step discrete dynamics
A = exp(A_SSP * T_SSP) * [1, T_DSP; 0, 1];
B = exp(A_SSP * T_SSP) * [-1; 0];

A = [2, -1; 
     1, 0.3];
B = [1; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION

% simulation time
t_max = 10;

% initial state
x0 = [0.5; 0.1];

% desired state
xd = [0.5; 0];

% MPC parameters
N = 5;
nx = size(A, 2);
nu = size(B, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZER

x_var = sdpvar(nx, N);    % state, [x0, x1, ..., xN]
u_var = sdpvar(nu, N-1);  % input, [u0, u1, ..., uN-1]

x0_ = sdpvar(nx, 1);      % initial state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINTS

% state constraint set, X = {x | H_x * x <= h_x}
x1_min = -1;
x1_max = 1;
x2_min = -1;
x2_max = 1;

H_x = [1,  0;
      -1,  0;
       0,  1;
       0, -1];
h_x = [x1_max;
      -x1_min;
       x2_max;
      -x2_min];

% input constraint set, U = {u | H_u * u <= h_u}
u_min = -0.1;
u_max = 0.1;

H_u = [1;
      -1];
h_u = [u_max;
      -u_min];

% initialize constraint list and add the constraint to it
x0_cons = x_var(:, 1) == x0_;
constraints = [x0_cons];

% add the state and input constraints at each node point
for i = 1 : N-1
    
    % new state constraint, x_i in Xcal
    state_cons = H_x * x_var(:, i) <= h_x;

    % new input constraint, u_i in Ucal
    input_cons = H_u * u_var(:, i) <= h_u;

    % new dynamics constraint
    dyns_cons = x_var(:, i+1) == A * x_var(:, i) + B * u_var(:, i);

    % add the constraints to the list
    constraints = [constraints, state_cons, input_cons, dyns_cons];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective function 

% state and input scaling
x_weights = [1; 1];
u_weights = [1];

% terminal state scaling
xN_weights = [1; 1];

% build the objective function
Qx = kron(eye(N), diag(x_weights));    % [Qx1, ..., QxN]
Qu = kron(eye(N-1), diag(u_weights));  % [Qu1, ..., QuN-1]
Q = blkdiag(Qx, Qu);                   % [Qx1, ..., QxN, Qu1, ..., QuN-1]
V = diag(xN_weights);                  % terminal cost

% build the objeciive function, stack x and u
e_var = x_var(:) - kron(ones(N, 1), xd);
objective = (1/2) * ([e_var; u_var(:)]' * Q * [e_var; u_var(:)]) + ...
            (1/2) * (x_var(:,N) - xd)' * V * (x_var(:,N) - xd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVER

% confgiure the optimizer settings
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% create the optimizer object
P = optimizer(constraints, objective, options, {x0_}, {x_var, u_var});

sol = P{x0}