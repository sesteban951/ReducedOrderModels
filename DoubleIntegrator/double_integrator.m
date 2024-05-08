%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics of the double intgerator system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% intial condition
x0 = [0; 0];

% time span
tspan = [0, 5];

% simulate the dynamics
[t, x] = ode45(@dynamics, tspan, x0);

% plot the state evolution
figure(1); hold on; grid on;
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2)
plot(x(1,1), x(1,2), 'rx', 'LineWidth', 2)
plot(x(end,1), x(end,2), 'gx', 'LineWidth', 2)
xlabel('$q$', 'Interpreter', 'latex')
ylabel('$\dot{q}$', 'Interpreter', 'latex')
legend('Trajectory', 'Initial Condition', 'Final Condition')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xdot = dynamics(t, x)

    % Define the system
    A = [0, 1;
         0, 0];
    B = [0;
         1];

    % control input
    u = 1;

    % dynamics
    xdot = A*x + B*u;
end