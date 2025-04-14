%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Double Integrator MPC Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% MPC parameters
mpc_params.dt = 0.05; % sampling time
mpc_params.N = 20;    % prediction horizon
mpc.u_lim = 5;        % double integrator control limit

% MPC costs
Q = diag([1, 1, 0.1, 0.1]);  % state cost
R = diag([0.01, 0.01]);      % input cost
V = diag([2, 2]);            % terminal cost
[H, f] = build_mpc_matrices(mpc_params, Q, R, V);

% initial state
x0 = [3;  % px
      2;  % py
      0;  % vx
      0]; % vy

% final state
xf = [0;  % px
      0;  % py
      0;  % vx
      0]; % vy

% simulation config
T_max = 10.0;
dt = mpc_params.dt;

% simulate the double integrator
t_size = T_max / mpc_params.dt;
T = zeros(t_size + 1, 1);
X = zeros(t_size + 1, 4);
U = zeros(t_size, 2);

% populate the intial condition
T(1) = 0;
X(1, :) = x0';

% simulate
t_increment = [0, dt];
xk = x0;
for i = 2:t_size+1

    % current time
    t = (i - 1) * dt;

    % compute the control input
    u = [cos(t); sin(t)];
    
    % update the state
    [t_flow, xk_flow] = ode45(@(t, x) dynamics(t, x, u), t_increment, xk);
    xk = xk_flow(end, :)';

    % store everything
    T(i) = t;
    X(i, :) = xk';
    U(i - 1, :) = u';
end

% plot thhe results
figure(1); 

subplot(2, 2, [1,3])
hold on; grid on;
plot(X(:, 1), X(:, 2), 'b', 'LineWidth', 2)

subplot(2, 2, 2)
hold on; grid on;
plot(T(1:end-1), U(:, 1), 'b', 'LineWidth', 2)
plot(T(1:end-1), U(:, 2), 'r', 'LineWidth', 2)

subplot(2, 2, 4)    
hold on; grid on;
plot(T(1:end-1), X(1:end-1, 3), 'b', 'LineWidth', 2)
plot(T(1:end-1), X(1:end-1, 4), 'r', 'LineWidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dybnamics of the double integrator system
function xdot = dynamics(~, x, u)
    
    % continuous parameters
    A = [zeros(2), eye(2);
         zeros(2), zeros(2)];
    B = [zeros(2); eye(2)];

    % dynamics
    xdot = A * x + B * u;
end

% build the MPC matrices
function build_mpc_matrices(mpc_params, Q, R, V)

    % 

end