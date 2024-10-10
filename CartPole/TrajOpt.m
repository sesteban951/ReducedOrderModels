clear all; close all; clc; 

% import casadi
import casadi.*

% Define the state and control variables
x = SX.sym('x',4); % state of the cart-pole
u = SX.sym('u');   % control input

% Define the continuous time dynamics
xdot = cartpole_dynamics(x, u);

% Define the objective function
L = objective_function(x, u);

% Continuous time dynamics function
f = Function('f', {x, u}, {xdot, L});

% Create an integrator function
T = 10;
N = 200;
dae = struct('x', x, 'p', u, 'ode', xdot, 'quad', L);
op = struct('tf', T/N, 'abstol', 1e-8); % integrator options
F = integrator('F', 'idas', dae, op);   % use the IDAS integrator

% Start with an empty NLP
w = {};
G = {};
J = 0;
lbw = [];
ubw = [];
lbg = [];
ubg = [];

% Initial conditions
Xk = MX.sym('X0', 4);
w{end+1} = Xk;
lbw = [lbw; 0.5; 0; 0; 0];  % initial state
ubw = [ubw; 0.5; 0; 0; 0];  % initial state

% Do it for the entire horizon
for k = 1:N
    
    % Local control
    Uk = MX.sym(['U' num2str(k-1)]);
    w{end+1} = Uk;
    lbw = [lbw; -200]; % lower bound on the control
    ubw = [ubw; 200];  % upper bound on the control

    % Call the integrator function 
    Fk = F('x0', Xk, 'p', Uk);
    Xk_end = Fk.xf;  % state at the end of the interval
    J = J + Fk.qf;   % accumulate cost

    % Local state at the next step
    Xk = MX.sym(['X' num2str(k)], 4);
    w{end+1} = Xk;
    lbw = [lbw; -inf; -inf; -inf; -inf]; % lower bounds on state
    ubw = [ubw; inf; inf; inf; inf];     % upper bounds on state

    % Dynamics constraint: Xk_end = Xk
    G{end+1} = Xk_end - Xk; 
end

% Create an NLP solver function
nlp = struct('x', vertcat(w{:}),...
             'f', J, 'g', vertcat(G{:}));
solver = nlpsol('solver', 'ipopt', nlp);

% Set solver options
opts = struct('print_time', 1, 'ipopt', struct('print_level', 0));
sol = solver('lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

% extract the state solution
x_opt = full(sol.x);
x1 = x_opt(1:5:end);
x2 = x_opt(2:5:end);
x3 = x_opt(3:5:end);
x4 = x_opt(4:5:end);
u = x_opt(5:5:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics taken from https://underactuated.mit.edu/acrobot.html#cart_pole

function L = objective_function(x, u)
    
    % desired states
    x1_des = 0;  % cart position
    x2_des = pi;  % pole angle
    x3_des = 0;  % cart velocity
    x4_des = 0;  % pole angular velocity
    
    % least squares cost
    Qx = diag([1, 30, 1, 1]);
    Qu = 1;

    % least sqaure objective function
    e = [x(1) - x1_des;
         x(2) - x2_des;
         x(3) - x3_des;
         x(4) - x4_des];
    L = e'*Qx*e + u'*Qu*u;
end

function xdot = cartpole_dynamics(x, u)

    % Constants
    mc = 1; % Mass of the cart
    mp = 0.2; % Mass of the pole
    g = 9.81; % Gravity
    l = 0.5; % Length of the pole

    % compute the accelerations
    xdot = [x(3);
            x(4);
            1/(mc + mp*sin(x(2))^2) * (u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2))));
            1/(l*(mc + mp*sin(x(2))^2)) * (-u*cos(x(2)) - mp*l*x(4)^2*cos(x(2))*sin(x(2)) - (mc + mp)*g*sin(x(2)))];
end