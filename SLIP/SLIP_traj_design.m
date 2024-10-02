%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP Traj Design (based on paper Wensing and Orin 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP paramss
params.m = 22;           % CoM mass (Achilles mass 22 kg)
params.g = 9.81;         % gravity
params.l0 = 0.5;         % spring free length (Achilles leg length 0.7 m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_running_params = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pick a desired forward velocity
v_des = 1.5;
t_stance_des = stance_time(v_des);

% TEST: query apex
x0 = [1.0, 0.0];    % pz_0, vx_0
u0 = [0.1, 15000];  % theta, ks
[t_stance, x_apex] = apex_return(x0, u0, params)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RUNNING PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_running_params == 1

    % choose range of desired velocities
    v_min = 0.0;
    v_max = 3.0;
    v = v_min: 0.01 :v_max;

    % plot the cadence
    subplot(2,1,1);
    c = cadence(v);
    plot(v, c, 'b', 'LineWidth', 2);
    xlabel('Desired Velocity [m/s]');
    ylabel('Cadence [steps/min]');
    title('Cadence vs. Desired Velocity');
    grid on;

    % plot the desired stance time
    subplot(2,1,2);
    t_s = stance_time(v);
    plot(v, t_s, 'r', 'LineWidth', 2);
    xlabel('Desired Velocity [m/s]');
    ylabel('Stance Time [s]');
    title('Stance Time vs. Desired Velocity');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cadence Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the cadence based on desired velocity, c [steps / min]
function c = cadence(v_des)
    c = 2.55 * v_des.^2 - 8.77 * v_des + 172.9;
end

% compute the stance time based on desired velocity
function t_stance = stance_time(v_des)
    % t_stance = -0.64 * log10(v_des) -0.2;
    t_stance = 10^(-0.2) * v_des.^(-0.82);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apex-to-Apex Dynamics Query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_stance, x_apex] = apex_return(x, u, params)

    % unpack the apex state
    z = x(1);
    vx = x(2);

    % unpack the control
    theta = u(1); % inpact angle
    ks = u(2);    % spring stiffness

    % sim params
    freq = 100;
    dt = 1/freq;
    tspan = 0:dt:3.0;  % to allow for switching before timeout

    % inital conditions
    x0 = [0;
          z;
          vx;
          0];
    domain = "flight";

    % assign the impact map
    alpha = theta;

    % set the switching manifolds
    options_f2g = odeset('Events', @(t,x)flight_to_ground(t, x, alpha, params), 'RelTol', 1e-8, 'AbsTol', 1e-9);
    options_g2f = odeset('Events', @(t,x)ground_to_flight(t, x, params), 'RelTol', 1e-8, 'AbsTol', 1e-9);
    options_poincare = odeset('Events', @(t,x)poincare_section(t, x), 'RelTol', 1e-8, 'AbsTol', 1e-9);
    
    % simulate the hybrid system
    t_current = 0;
    num_transitions = 0;
    max_num_transitions = 2;
    T_apex = [];  % apex time storage
    X_apex = [];  % apex state storage
    T_TD = [];  % touch down times
    T_LO = [];  % lift off times

    % forward propagate the sim
    while num_transitions <= max_num_transitions
        
        % switch domains
        if domain == "flight"
            
            % flight: x = [x, z, x_dot, z_dot]
            [t_flight, x_flight] = ode45(@(t,x)dynamics_f(t,x,params), tspan + t_current, x0, options_f2g);
            [~, ~, t_apex, x_apex, ~] = ode45(@(t,x)dynamics_f(t,x,params), tspan + t_current, x0, options_poincare); % purely used for poincare section
            T_apex = [T_apex; t_apex];
            X_apex = [X_apex; x_apex];

            % store the trajectory
            T_TD = [T_TD; t_flight(end)];
    
            % udpate the current time and the intial state
            t_current = t_flight(end);

            % compute the foot position
            p_foot = [x_flight(end,1) + params.l0 * sin(alpha); 
                      x_flight(end,2) - params.l0 * cos(alpha)];

            % set new initial condition
            x0 = x_flight(end,:);
            x0 = cart_to_polar(x0, params, alpha);

            % define new domain
            domain = "ground";
            num_transitions = num_transitions + 1;

        elseif domain == "ground"
            
            % ground: x = [r, theta, r_dot, theta_dot]
            [t_ground, x_ground] = ode45(@(t,x)dynamics_g(t,x,params, ks), tspan + t_current, x0, options_g2f); 

            % convert the polar state to cartesian
            for i = 1:length(t_ground)
                x_ground(i,:) = polar_to_cart(x_ground(i,:)); % convert it to cartesian
                x_ground(i,1) = x_ground(i,1) + p_foot(1);    % add the foot position offset
            end

            % store the trajectory
            T_LO = [T_LO; t_ground(end)];

            % udpate the current time and the intial state
            t_current = t_ground(end);

            % set new initial condition
            x0 = x_ground(end,:);
        
            % define new domain
            domain = "flight";
            num_transitions = num_transitions + 1;
        end
    end

    t_stance = T_TD(end) - T_LO(end);
    x_apex = X_apex(end,:);

end

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
function xdot = dynamics_g(t, x_polar, params, k)
    
    % unpack the system parameters
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

% guard: flight to ground
function [value, isterminal, direction] = flight_to_ground(t_abs, x_cart, alpha, params)

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

% apex: detect the peak of the flight phase
function [value, isterminal, direction] = poincare_section(~, x_cart)

    % poincare section
    z_dot = x_cart(4);

    value = z_dot;  % z_dot = 0
    isterminal = 0; % stop integrating
    direction = -1; % negative direction

end

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