%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spinrg Loaded Inverted Pendulum (SLIP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP paramss
params.k = 250;  % spring constant [N/m]
params.m = 1;    % CoM mass (Achilles mass 22 kg)
params.g = 9.81;  % gravity
params.l0 = 0.6;  % spring free length (Achilles leg length 0.7 m)
params.K = 0.20;  % Raibert controller gain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sim params
dt = 0.01;
tspan = 0:dt:3.0;  % to allow for switching before timeout

% initial conditions (always start in flight)
x0 = [0.0;   % x
      1.0;   % z
      0.1;  % x_dot
      0.0];  % z_dot
domain = "flight";

% initial foot angle
alpha_prev = 0;
alpha = angle_control(x0, params);

% set the switching manifolds
options_f2g = odeset('Events', @(t,x)flight_to_ground(t, x, params), 'RelTol', 1e-7, 'AbsTol', 1e-8);
options_g2f = odeset('Events', @(t,x)ground_to_flight(t, x, params), 'RelTol', 1e-7, 'AbsTol', 1e-8);

% simulate the hybrid system
t_current = 0;
num_transitions = 0;
max_num_transitions = 25;
D = [];  % domain storage
T = [];  % time storage
X = [];  % state storage
F = [];  % ground foot position storage
while num_transitions <= max_num_transitions
    
    % switch domains
    if domain == "flight"
        
        disp("flight")

        % flight: x = [x, z, x_dot, z_dot]
        [t_flight, x_flight] = ode45(@(t,x)dynamics_f(t,x,params), tspan, x0, options_f2g);

        % store the trajectory
        D = [D; 0 * ones(length(t_flight),1)];
        T = [T; t_flight + t_current];
        X = [X; x_flight];
  
        % udpate the current time and the intial state
        t_current = T(end);

        % compute foot trajectory
        for i = 1:length(t_flight)
            p_foot = [x_flight(i,1) + params.l0 * sin(alpha); 
                      x_flight(i,2) - params.l0 * cos(alpha)];
            F = [F; p_foot'];
        end

        % set new initial condition
        x0 = x_flight(end,:);
        x0 = cart_to_polar(x0, params, alpha);

        % define new domain
        domain = "ground";
        num_transitions = num_transitions + 1;

    elseif domain == "ground"
        
        disp("ground")

        % ground: x = [r, theta, r_dot, theta_dot]
        [t_ground, x_ground] = ode45(@(t,x)dynamics_g(t,x,params), tspan, x0, options_g2f); 

        % convert the polar state to cartesian
        for i = 1:length(t_ground)
            x_ground(i,:) = polar_to_cart(x_ground(i,:)); % convert it to cartesian
            x_ground(i,1) = x_ground(i,1) + p_foot(1);            % add the foot position offset
        end

        % store the trajectory
        D = [D; 1 * ones(length(t_ground),1)];
        T = [T; t_ground + t_current];
        X = [X; x_ground];

        % udpate the current time and the intial state
        t_current = T(end);

        % compute foot trajectory
        for i = 1:length(t_ground)
            p_foot = F(end,:);
            F = [F; p_foot];
        end

        % set new initial condition
        x0 = x_ground(end,:);
        alpha = angle_control(x0, params);
    
        % define new domain
        domain = "flight";
        num_transitions = num_transitions + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a new figure
figure('Name', 'SLIP Simulation');
hold on;
yline(0);
xlabel('$p_x$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$p_z$ [m]', 'Interpreter', 'latex', 'FontSize', 16);

% set axis limits
x_min = min(X(:,1)) - 0.1;
x_max = max(X(:,1)) + 0.1;
z_min = -0.1;
z_max = max(X(:,2)) + 0.1;
xlim([x_min, x_max]);
ylim([z_min, z_max]);

tic;
t_now = T(1);
ind = 1;
while t_now < T(end)

    % plot the foot and pole
    pole = plot([F(ind,1), X(ind,1)], [F(ind,2), X(ind,2)], 'k', 'LineWidth', 2.5);
    foot = plot(F(ind,1), F(ind,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');  % in flight

    % plot the SLIP COM
    if D(ind) == 0
        com = plot(X(ind,1), X(ind,2), 'ko', 'MarkerSize', 20, 'MarkerFaceColor', 'b');  % on the ground
    elseif D(ind) == 1
        com = plot(X(ind,1), X(ind,2), 'ko', 'MarkerSize', 20, 'MarkerFaceColor', 'r');  % in flight
    end
    
    % current time
    time = sprintf("Time = %.2f", T(ind));
    title(time,'Interpreter','latex', 'FontSize', 16)   
    
    drawnow;

    % wait until the next time step
    while toc < T(ind+1)
        % wait
    end

    % increment the index
    if ind+1 == length(T)
        break
    else
        ind = ind + 1;
        delete(pole);
        delete(foot);
        delete(com);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SLIP flight dynamics
function xdot = dynamics_f(~, x_cart, params)
    
    % cartesian state, x = [x, z, x_dot, z_dot]
    x_dot = x_cart(3);
    z_dot = x_cart(4);

    % drift dynamics
    xdot = [x_dot;
            z_dot;
            0;
            -params.g];
end

% SLIP ground dynamics
function xdot = dynamics_g(~, x_polar, params)
    
    % unpack the system parameters
    k = params.k;
    m = params.m;
    g = params.g;
    l0 = params.l0;

    % polar state, x = [r, theta, r_dot, theta_dot]
    r = x_polar(1);
    theta = x_polar(2);
    r_dot = x_polar(3);
    theta_dot = x_polar(4);

    xdot = [r_dot;
            theta_dot;
            r*theta_dot^2 - g*cos(theta) + (k/m)*(l0 - r);
            -(2/r)*r_dot*theta_dot + (g/r)*sin(theta)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control input
function alpha = angle_control(x_cart, params)
    
    % simple Raibert controller
    K = params.K;
    alpha = K * x_cart(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUARDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guard: flight to ground
function [value, isterminal, direction] = flight_to_ground(~, x_cart, params)

    % get the angle of the foot
    alpha = angle_control(x_cart, params);

    % to determine if the SLIP foot has hit the ground
    z_com = x_cart(2);
    foot_height = z_com - params.l0 * cos(alpha);

    % guard conditions
    value = foot_height;  % foot height at ground
    isterminal = 1;       % 1: stop integrating
    direction = -1;       % direction 
end

% guard: ground to flight, leg condition
function [value, isterminal, direction] = ground_to_flight(~, x_polar, params)

    % equivalent representation in cartesian coordinates
    x_cart = polar_to_cart(x_polar);

    % leg length is uncompressed, r = l0
    l = [x_cart(1); x_cart(2)];
    leg_length = norm(l);                        % Euclidean length of the leg
    compressed_length = leg_length - params.l0;  % difference from nominal uncompressed length

    % taking off condition, vel >= 0
    xdot = x_cart(3); 
    zdot = x_cart(4);  
    v_com = [xdot; zdot];  
    vel = l' * v_com;        % velocity of the CoM along the leg direction

    if compressed_length >= 0
        if vel >= 0
            value = 0;
        else
            value = 1;
        end
    else
        value = 1;
    end

    % Ensure the solver stops when both conditions are met
    isterminal = 1;  % Stop integration when event is triggered
    direction =  0;  % compressed_length must be increasing (zero-crossing from positive to negative)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert caterisan <---> polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert a cartesian state to polar state where the origin is at the foot
function x_polar = cart_to_polar(x_cart, params, alpha)
    
    % flight state, x = [x, z, x_dot, z_dot]
    x = x_cart(1);
    z = x_cart(2);
    xdot = x_cart(3);
    zdot = x_cart(4);

    % positions
    p_com = [x; z];  % CoM position
    p_foot = [x + params.l0 * sin(alpha); z - params.l0 * cos(alpha)]; % foot position

    x = p_com(1) - p_foot(1);
    z = p_com(2) - p_foot(2);

    r = sqrt(x^2 + z^2);
    th = atan2(x, z);     % be carefule about arctan2
    rdot = (x*xdot + z*zdot) / r;
    thdot = (xdot*z - x*zdot) / r^2;

    x_polar = [r; th; rdot; thdot];
end

% convert a polar state to cartesian state, where the origin is at the foot
function x_cart = polar_to_cart(x_polar)
    
    % ground state, x = [r, theta, r_dot, theta_dot]
    r = x_polar(1);
    th = x_polar(2);
    rdot = x_polar(3);
    thdot = x_polar(4);

    x = r * sin(th);
    z =  r * cos(th);
    xdot = (rdot * sin(th) + r * thdot * cos(th));
    zdot =  rdot * cos(th) - r * thdot * sin(th);

    x_cart = [x; z; xdot; zdot];
end
