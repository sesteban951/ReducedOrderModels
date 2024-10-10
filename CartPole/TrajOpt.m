%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CART POLE Multiple-Shooting Trajectory Optimization w/ CasADi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; 

% import casadi
import casadi.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System Parameters
params.mc = 1;   % Mass of the cart
params.mp = 0.2; % Mass of the pole
params.g = 9.81; % Gravity
params.l = 0.5;  % Length of the pole

% instatiate the solver parameters
T = 15;  % total time horizon
N = 300; % total number of nodes

% quadratic weights
Q = diag([1, 30, 1, 1]); % state weights
R = 1;                   % control weights

% bounds on the state for all time steps
x0 = [0;  % cart pos
      0;  % pole angle
      0;  % cart vel
      0]; % pole angular vel
xd = [0;  % cart pos
      pi;  % pole angle
      0;  % cart vel
      0]; % pole angular vel
x_min = [-inf; -inf; -inf; -inf];  
x_max = [ inf;  inf;  inf;  inf];
u_min = -200;
u_max =  200;

% overall options
dyn_integrator = 'idas'; % integrator type
nlp_solver = 'ipopt';   % nlp solver type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the state and control variables
x = SX.sym('x',4); % state of the cart-pole
u = SX.sym('u');   % control input

% Define the continuous time dynamics
xdot = cartpole_dynamics(x, u, params);

% Define the objective function
L = objective_function(x, u, xd, Q, R);

% Create an integrator function
dyn = struct('x', x, 'p', u, 'ode', xdot, 'quad', L); % instantiate the dynamics
op = struct('tf', T/N, 'abstol', 1e-8);               % integrator options
F = integrator('F', dyn_integrator, dyn, op);         % use the IDAS integrator

% Start with empty NLP constraints
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
lbw = [lbw; x0];  % initial state
ubw = [ubw; x0];  % initial state

% Do it for the entire horizon
for k = 1:N
    
    % Local control
    Uk = MX.sym(['U' num2str(k-1)]);
    w{end+1} = Uk;
    lbw = [lbw; u_min]; % lower bound on the control
    ubw = [ubw; u_max];  % upper bound on the control

    % Call the integrator function 
    Fk = F('x0', Xk, 'p', Uk);
    Xk_end = Fk.xf;  % state at the end of the interval
    J = J + Fk.qf;   % accumulate cost

    % Local state at the next step
    Xk = MX.sym(['X' num2str(k)], 4);
    w{end+1} = Xk;
    lbw = [lbw; x_min]; % lower bounds on state
    ubw = [ubw; x_max];     % upper bounds on state

    % Dynamics constraint: Xk_end = Xk
    G{end+1} = Xk_end - Xk; 
end

% Create an NLP solver function
nlp = struct('x', vertcat(w{:}),...
             'f', J, 'g', vertcat(G{:}));
solver = nlpsol('solver', nlp_solver, nlp);

% Set solver options
opts = struct('print_time', 1, nlp_solver, struct('print_level', 0));
sol = solver('lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

% extract the state solution
x_opt = full(sol.x);
p = x_opt(1:5:end);      % cart position
th = x_opt(2:5:end);     % pole angle
pdot = x_opt(3:5:end);  % cart velocity
thdot = x_opt(4:5:end); % pole angular velocity
u = x_opt(5:5:end);      % control input

% save the data into csv files
X = [p, th, pdot, thdot];
U = u;
T = linspace(0, T, N+1)';

csvwrite('data/state_data.csv', X);
csvwrite('data/input_data.csv', U);
csvwrite('data/time_data.csv', T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% least squares loss
function L = objective_function(x, u, xd, Q, R)

    % least sqaure objective function
    e = [x(1) - xd(1);
         x(2) - xd(2);
         x(3) - xd(3);
         x(4) - xd(4)];
    L = e' * Q * e + u' * R * u;
end

% Dynamics taken from https://underactuated.mit.edu/acrobot.html#cart_pole
function xdot = cartpole_dynamics(x, u, params)

    % system parameters
    mc = params.mc;
    mp = params.mp;
    g = params.g;
    l = params.l;

    % the dynamics
    xdot = [x(3);
            x(4);
            1/(mc + mp*sin(x(2))^2) * (u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2))));
            1/(l*(mc + mp*sin(x(2))^2)) * (-u*cos(x(2)) - mp*l*x(4)^2*cos(x(2))*sin(x(2)) - (mc + mp)*g*sin(x(2)))];
end