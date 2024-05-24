%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Reachbility of the inverted pendulum example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; 

% System dynamics
params.m = 1.0;  % pendulum mass
params.l = 0.5;  % pendulum length
params.b = 0.2;  % pendulum damping
params.g = 9.81; % gravity

% example ODE solver
x0 = [pi/4;  % pole position
      0.01];    % pole velocity

% get linearization functions
[f, g, Df, Dg] = get_functions(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation time
dt = 0.01;
t_max = 5;
tspan = 0 : dt : t_max;

% linearize at the intial condition and solve the ODE
[A, B, C] = get_linearization(x0, 0, f, g, Df, Dg);
x_t = solve_LTI(A, B, C, tspan, x0);

% solve the ODE
[t, x] = ode45(@(t,x) dynamics(t,x, params), tspan, x0);

% plots in state space
figure(); hold on; grid on;

% find the limits of the state space
x1_min = min(x(:,1));
x1_max = max(x(:,1));
x2_min = min(x(:,2));
x2_max = max(x(:,2));
xlims = [x1_min, x1_max] * 1.1;
ylims = [x2_min, x2_max] * 1.1;
[Theta, Theta_dot] = meshgrid(linspace(xlims(1), xlims(2), 50), linspace(ylims(1), ylims(2), 50));

% plot the vector field
U = zeros(size(Theta));
V = zeros(size(Theta_dot));
for i = 1:numel(Theta)
    xdot = dynamics(0, [Theta(i); Theta_dot(i)], params);
    U(i) = xdot(1);
    V(i) = xdot(2);
end
quiver(Theta, Theta_dot, U, V, 'k');

% plot the state space
plot(x(:,1), x(:,2), 'b', 'LineWidth', 2); 
plot(x_t(:,1), x_t(:,2), 'r', 'LineWidth', 2);
xlim(xlims);
ylim(ylims);
xline(0);
yline(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% noninear dynamics
function xdot = dynamics(t, x, params)

    % state
    theta = x(1);
    theta_dot = x(2);
    
    % control input
    u = 0;
    
    % dynamics
    f = [theta_dot;
        -params.g/params.l * sin(theta) - params.b/(params.m*params.l^2) * theta_dot];
    g = [0; 
         1/(params.m*params.l^2)];
    xdot  = f + g * u;
end

% solve the linear dynamics
function x = solve_LTI(A, B, C, tspan, x0)

    % solve the ODE
    x = zeros(length(tspan), 2);
    for i =1:length(tspan)
        x_t = expm(A*tspan(i)) * x0 + (expm(A*tspan(i)) - eye(size(A,2))) * inv(A) * C;  
        x(i,:) = x_t';
    end
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

    
