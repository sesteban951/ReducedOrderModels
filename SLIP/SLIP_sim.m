%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spinrg Loaded Inverted Pendulum (SLIP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP params
params.k = 7500; % spring constant [N/m]
params.m = 22;    % CoM mass (Achilles mass 22 kg)
params.g = 9.81;  % gravity
params.l0 = 0.6;  % spring free length (Achilles leg length 0.7 m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sim params
dt = 0.005;
tspan = 0:dt:2.0;  % to allow for switching before timeout

% initial conditions (always start in flight)
x0 = [0.0;   % x
      1.0;   % z
      0.0;   % x_dot
      0.0];  % z_dot
domain = "flight";

% set the switching manifolds
options_f2g = odeset('Events', @(t,x)flight_to_ground(t, x, params), 'RelTol', 1e-5, 'AbsTol', 1e-7);
options_g2f = odeset('Events', @(t,x)ground_to_flight(t, x, params), 'RelTol', 1e-5, 'AbsTol', 1e-7);

% simulate the hybrid system
t_current = tspan(1);
num_transitions = 0;
max_num_transitions = 2;
T = [];
X = [];
while num_transitions <= max_num_transitions
    
    % switch domains
    if domain == "flight"
        
        disp("flight")
        
        % flight: x = [x, z, x_dot, z_dot]
        [t_flight, x_flight] = ode45(@(t,x)dynamics_f(t,x,params), tspan, x0, options_f2g);

        % store the trajectory
        T = [T; t_flight + t_current];
        X = [X; x_flight];
  
        figure; hold on;
        plot(t_flight, x_flight(:,2), 'b.');
        xlabel('x'); ylabel('z');
        yline(0);

        % udpate the current time and the intial state
        t_current = T(end);
    
        % apply reset map (includes cartesian to polar conversion)
        x0 = x_flight(end,:);
        x0 = cart_to_polar(x0)

        % define new domain
        domain = "ground";
        num_transitions = num_transitions + 1;

    elseif domain == "ground"
        
        disp("ground")

        % ground: x = [r, theta, r_dot, theta_dot]
        [t_ground, x_ground] = ode45(@(t,x)dynamics_g(t,x,params), tspan, x0, options_g2f); 

        % store the trajectory
        T = [T; t_ground + t_current];
        X = [X; x_ground];

        figure; hold on;
        x_pos = -x_ground(:,1) .* sin(x_ground(:,2));
        z_pos = x_ground(:,1) .* cos(x_ground(:,2));
        plot(t_ground, z_pos, 'r.');

        % udpate the current time and the intial state
        t_current = T(end);

        % apply reset map (includes polar to cartesian conversion)
        x0 = x_ground(end,:);
        x0 = polar_to_cart(x0)
    
        % define new domain
        domain = "flight";
        num_transitions = num_transitions + 1;
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
% GUARDS and RESET MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guard: flight to ground
function [value, isterminal, direction] = flight_to_ground(~, x_cart, params)
    
    % to determine if the SLIP foot has hit the ground
    z_com = x_cart(2);
    l0 = params.l0;
    foot_height = z_com - l0;

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
    leg_length = norm(l);  % Euclidean length of the leg
    compressed_length = leg_length - params.l0;  % difference from nominal uncompressed length

    % taking off condition, vel >= 0
    xdot = x_cart(3); 
    zdot = x_cart(4);  
    v_com = [xdot; zdot];  
    l_unit = l / leg_length;       % unit vector along the leg
    v_unit = v_com / norm(v_com);  % unit vector of the CoM velocity
    vel = l_unit' * v_unit;        % velocity of the CoM along the leg direction

    if compressed_length >= 0
        if vel >= 0
            value = 0;
        end
    else
        value = compressed_length;
    end

    % Ensure the solver stops when both conditions are met
    isterminal = 1;  % Stop integration when event is triggered
    direction =  0;  % compressed_length must be decreasing (zero-crossing from positive to negative)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert caterisan <---> polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert a cartesian state to polar state where the origin is at the foot
function x_polar = cart_to_polar(x_cart)
    x = x_cart(1);
    z = x_cart(2);
    xdot = x_cart(3);
    zdot = x_cart(4);

    r = sqrt(x^2 + z^2);
    th = atan2(x, z);     % be carefule about arctan2
    rdot = (x*xdot + z*zdot) / r;
    thdot = (xdot*z - x*zdot) / r^2;

    x_polar = [r; th; rdot; thdot];
end

% convert a polar state to cartesian state, where the origin is at the foot
function x_cart = polar_to_cart(x_polar)
    r = x_polar(1);
    th = x_polar(2);
    rdot = x_polar(3);
    thdot = x_polar(4);

    y = -r * sin(th);
    z =  r * cos(th);
    ydot = -rdot * sin(th) - r * thdot * cos(th);
    zdot =  rdot * cos(th) - r * thdot * sin(th);

    x_cart = [y; z; ydot; zdot];
end
