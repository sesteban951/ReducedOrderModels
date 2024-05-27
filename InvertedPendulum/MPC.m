%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive MPC controller for the inverted pendulum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% System dynamics
params.m = 1.0;  % pendulum mass
params.l = 1.0;  % pendulum length
params.b = 0.1;  % pendulum damping
params.g = 9.81; % gravity

% get linearization functions
[f, g, Df, Dg] = get_functions(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intial and desired states
x0 = [0;    % pole position
      0];    % pole velocity

xd = [pi;
      0]; % desired state

% total sim time
t_span = 7;

% specify the MPC parameters
dt = 0.04;
N = 80;

% specify the state and input vector sizes
nx = 2;
nu = 1;

% state and input size for MPC horizon
nx_horizon = nx * (N-1);
nu_horizon = nu * (N-1);

% add terminal constraint
terminal_constraint = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear the optimization variables
yalmip('clear');               % good practice to clear the variables

% state and input variables
states = sdpvar(nx, N);     % state variables (nx x N)
inputs = sdpvar(nu, N-1);   % input variables (nu x N-1)

x_init = sdpvar(nx, 1);     % initial state variable (nx x 1)

% containers to hold all linearized discrete dynamics matrices
Ad = sdpvar(nx, nx_horizon);  % linear discrete drift matrix
Bd = sdpvar(nx, nu_horizon);  % linear discrete actuation matrix
Cd = sdpvar(nx, N-1);         % linear discrete constant term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the constraints list
% state constraints, x in X
th_max =  5;          % pole position constraints
th_min = -5;
th_dot_max =  30;       % pole velocity conastraints
th_dot_min = -30;

% input constraints, u in U
u_max =  20;
u_min = -20;

% build the conastriants as a H-polyhedron, C(x, u) <= [d_x; d_u]
C_x = kron(eye(nx), [1; -1]);
d_x = [th_max; 
      -th_min; 
       th_dot_max; 
      -th_dot_min];
    
C_u  = kron(eye(nu), [1; -1]);
d_u = [u_max; 
      -u_min];

% initialize the constraint list and add the state and input constraints
init_const = (states(:,1) == x_init);
constraints = [init_const];

% for loop to add all the constraints
for i = 1 : N-1
    
    % new state constraint
    state_cons = (C_x * states(:,i) <= d_x);

    % new input constraint
    input_cons = (C_u * inputs(:,i) <= d_u);

    % new dynamics constraint
    dyn_cons = states(:, i+1) ==  Ad(:, (i-1)*nx + 1 : i*nx) * states(:, i) ...
                                + Bd(:,i) * inputs(:,i) ...
                                + Cd(:,i);

    % add the constraints
    constraints = [constraints; state_cons; input_cons; dyn_cons];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the objective function

% weights for the objective function
x_weights = [5;          % pole posititon
             1];         % pole velocity
x_weights_terminal = [5; % pole posititon
                       1];  % pole velocity
u_weights = [0.05];       % control effort

% build Qx and Qu matrices and add to make big Q matrix
Qx = diag(x_weights);
Qu = diag(u_weights);
x_cost = kron(eye(N), Qx);
u_cost = kron(eye(N-1), Qu);
Q = blkdiag(x_cost, u_cost);

% build the objective function with terminal cost, but add no terminal constraint
xd_horizon = kron(ones(N, 1), xd);
if (terminal_constraint == true)

    % build the objective function with terminal constraint
    objective = (1/2) * ([states(:)-xd_horizon; inputs(:)]' * Q * [states(:)-xd_horizon; inputs(:)]);

    %  add the terminal constraint
    term_cons = (states(:,N) == xd);
    constraints = [constraints; term_cons];
else

    % build the objective function with terminal cost
    V =  diag(x_weights_terminal);
    objective = (1/2) * ([states(:)-xd_horizon; inputs(:)]' * Q * [states(:)-xd_horizon; inputs(:)]) + ...
                (1/2) * (states(:,N)-xd)' * V * (states(:,N)-xd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the optimizer

% options for the solver
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% create the optimzer/function
P = optimizer(constraints, ...             % constraints
              objective, ...               % objective
              options, ...                 % options
              {Ad, Bd, Cd, x_init}, ...    % input parameters
              {states, inputs});           % decision varaibles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC Time, Letssss GOOOO!!!

% to store the actual states and input trajectories
T = [0];    % history of time steps
X = [x0'];  % history of states
U = [];     % history of inputs

% forward simulate
u_bar = zeros(nu, N-1);
for t = 0:dt:t_span

    disp(['Time: ', num2str(t)]);  

    % forward simulate the nonlinar dynamics for the horizon (w/ some nominal control)
    % to get a nominal trajectory
    [t_rollout, x_bar] = dynamics_w_input(dt, x0, u_bar, params);
    x_bar = x_bar';

    % linearize the dynamics about this nominal trajectory
    Ad_list = [];
    Bd_list = [];
    Cd_list = [];

    % compute the linearization about this nomincal trajectory (x_nom, u_nom)
    for k = 1:N-1
        % nominal trajectory
        x_k = x_bar(:,k);
        u_k = u_bar(:,k);

        % compute the linearization for linear contiuous dynamics
        [Ac_k, Bc_k, Cc_k] = get_linearization(x_k, u_k, f, g, Df, Dg);

        % compute the discrete version of the linearized dynamics
        [Ad_k, Bd_k, Cd_k] = css2dss('Exact', dt, Ac_k, Bc_k, Cc_k);

        % store the linearization
        Ad_list = [Ad_list, Ad_k];
        Bd_list = [Bd_list, Bd_k];
        Cd_list = [Cd_list, Cd_k];
    end

    % solve the MPC problem for x_star and u_star given the linearized dynamics and x0
    [sol, diagnostics] = P({Ad_list, Bd_list, Cd_list, x0});  
    
    % extract the solution
    x_opt = sol{1};  % optimal state trajectory
    u_opt = sol{2};  % optimal input trajectory

    % take the first input and use it as a feedforward control
    u_ff = u_opt(:,1);

    % forward simulate the nonlinear dynamics with the feedforward control
    [~, x] = ode45(@(t, x) dynamics_dt(t, x, u_ff, params), [0:dt:dt], x0);

    % check if any NaN values
    if any(isnan(x0))
        disp('NaN values in the state trajectory. Stopping the simulation at time: ');
        disp(t);
        break;
    else

        % store the actual states and input trajectories
        T = [T, t];
        X = [X; x(end,:)];
        U = [U; u_ff'];

        % update the initial state
        x0 = x(end,:)';
    end
end

% animate using the results
figure(1); axis equal;
xlim = [-params.l-0.5, params.l+0.5];
ylim = [-params.l-0.5, params.l+0.5];
axis([xlim, ylim]);
pause(0.5);
tic;
for i = 1:length(T)-1
    while toc < T(i)
        % do nothing
    end
    
    % show the title
    title(sprintf('Time: %0.2f', T(i)));

    % draw the pendulum
    draw_pendulum(T(i), X(i,:), params);

    % clear the plot
    if i ~= length(T)-1
        delete(findobj('type','line'));
        delete(findobj('type','rectangle'));
    end
end

% plot the results
pause(1.0);
figure();

subplot(2,1,1); hold on;
plot(T, X(:,1), 'b', 'LineWidth', 2);
plot(T, X(:,2), 'r', 'LineWidth', 2);
xlabel('time');
ylabel('states');
yline(xd(1), '--', 'LineWidth', 2);
yline(xd(2), '--', 'LineWidth', 2);
yline(0);
legend('theta', 'theta dot');

subplot(2,1,2); hold on;
stairs(T(1:end-1), U, 'LineWidth', 2);
xlabel('time');
ylabel('input');
yline(u_max, '--', 'LineWidth', 2);
yline(u_min, '--', 'LineWidth', 2);
yline(0)
legend('u');

% plot the state trajectory
figure(); hold on;
plot(X(:,1), X(:,2), 'k', 'LineWidth', 2);
plot(xd(1), xd(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
xline(0);
yline(0);
xlabel('theta');
ylabel('theta dot');

disp("Done.")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nonlinear dynamics with inputs
function [T, X] = dynamics_w_input(dt, x0, u_list, params)

    % number of inputs 
    u_len = size(u_list, 2);

    % containers to hold the solution 
    X = zeros(u_len+1, length(x0));
    X(1, :) = x0';
    T = 0:dt:dt*u_len;

    % propagate the dynamics with zero order holds on the inputs
    for i = 1:length(u_list)
        
        % get the zero order hold to apply
        u = u_list(:,i);

        % apply zero order hold with the given initial condition
        [~, x_t] = ode45(@(t, x) dynamics_dt(t, x0, u, params), 0:dt:dt, x0);
        
        % save the solution
        x = x_t(end, :)';
        X(i+1, :) = x;
        
        % update the initial condition
        x0 = x';
    end

end

% poropagate the dynamics for dt and one zero order hold u
function x_dot = dynamics_dt(t, x, u, params)
    % dynamics
    theta = x(1);
    theta_dot = x(2);

    % drift and actuation matrices
    f_x = [theta_dot;
          -params.g/params.l * sin(theta) - params.b/(params.m*params.l^2) * theta_dot];
    g_x = [0; 
         1/(params.m*params.l^2)];

    % compute the dynamics
    x_dot = f_x + g_x * u;
end

% get the linearized dynamics at a point (x_bar, u_bar)
function [A, B, C] = get_linearization(x_bar, u_bar, f, g, Df, Dg)

    % compute A
    f_x = f(x_bar(1), x_bar(2));
    g_x = g;
    df_dx_x = Df(x_bar(1), x_bar(2));
    dg_dx_x = Dg(x_bar(1), x_bar(2));

    % compute A
    A = df_dx_x + dg_dx_x * u_bar;

    % compute B
    B = g_x;

    % compute C
    C = (f_x + g_x * u_bar) - (df_dx_x + dg_dx_x * u_bar) * x_bar - g_x * u_bar;
end

% get nonlin dynamics functions and linearizations of dynamics function, Df, Dg, for a point x
function [f, g, Df, Dg] = get_functions(params)

    % symbolic variables
    syms theta theta_dot

    % dynamics
    f_x = [theta_dot;
          -params.g/params.l * sin(theta) - params.b/(params.m*params.l^2) * theta_dot];
    g_x = [0; 
         1/(params.m*params.l^2)];

    % convert the drift and actuation vectors as symbolic functions
    f =  matlabFunction(f_x, 'Vars', {'theta', 'theta_dot'});
    g =  g_x;

    % compute the jacobians
    Df = [diff(f_x, theta), diff(f_x, theta_dot)];
    Dg = [diff(g_x, theta), diff(g_x, theta_dot)];

    % convert to function handles
    Df = matlabFunction(Df, 'Vars', {'theta', 'theta_dot'});
    Dg = matlabFunction(Dg, 'Vars', {'theta', 'theta_dot'});
end

    
