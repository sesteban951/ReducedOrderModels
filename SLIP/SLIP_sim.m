%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spinrg Loaded Inverted Pendulum (SLIP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP params
params.k = 150;   % spring constant
params.m = 22;    % CoM mass (Achilles mass)
params.g = 9.81;  % gravity
params.l0 = 0.7;  % spring free length (Achilles leg length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sim params
dt = 0.01;
tspan = 0:dt:100.0;  % to allow for switching before timeout

% initial conditions (always start in flight)
x0 = [0.0;   % x
      2.0;   % z
      0.5;   % x_dot
      0.0];  % z_dot
domain = "flight";

% set the switching manifolds
options_f2g = odeset('Events', @(t,x)flight_to_ground(t, x, params));

% simulate the hybrid system
t_current = tspan(1);
num_transitions = 0;
max_num_transitions = 8;
T = [];
X = [];
while num_transitions < max_num_transitions
    
    % switch domains
    if domain == "flight"
        
        disp("flight")
        
        % flight: x = [x, z, x_dot, z_dot]
        [t_flight, x_flight] = ode45(@(t,x)dynamics_f(t,x,params), tspan, x0, options_f2g);

        % store the trajectory
        T = [T; t_flight + t_current];
        X = [X; x_flight];

        % udpate the current time and the intial state
        t_current = T(end)
    
        % apply reset map
        x0 = x_flight(end,:);
        x0(4) = -0.8*x0(4); % flip the z_dot

        % define new domain
        domain = "ground";
        num_transitions = num_transitions + 1;

    elseif domain == "ground"
        
        disp("ground")

        % ground: x = [r, theta, r_dot, theta_dot]
        [t_ground, x_ground] = ode45(@(t,x)dynamics_f(t,x,params), tspan, x0, options_f2g);

        % store the trajectory
        T = [T; t_ground + t_current];
        X = [X; x_ground];

        % udpate the current time and the intial state
        t_current = T(end)

        % apply reset map
        x0 = x_ground(end,:);
        x0(4) = -0.8*x0(4); % flip the theta_dot

        % define new domain
        domain = "flight";
        num_transitions = num_transitions + 1;
    end

end

% plot the trajectory
figure;
plot(X(:,1), X(:,2)); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SLIP flight dynamics
function xdot = dynamics_f(~, x_cart, params)
    
    % cartesian state, x = [x, z, x_dot, z_dot]
    % x = x_cart(1);
    % z = x_cart(2);
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
function [value, isterminal, direction] = flight_to_ground(t, x_cart, params)
    
    % to determine if the SLIP foot has hit the ground
    z_com = x_cart(2);
    l0 = params.l0;
    foot_height = z_com - l0;

    % guard conditions
    value = foot_height;      % foot height at ground
    isterminal = 1; 
    direction = -1; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert caterisan <---> polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert a cartesian state to polar state where the origin is at the foot
% https://math.stackexchange.com/questions/2444965/relationship-between-cartesian-velocity-and-polar-velocity
function x_polar = cart_to_polar(x_cart)
    y = x_cart(1);
    z = x_cart(2);
    ydot = x_cart(3);
    zdot = x_cart(4);

    r = sqrt(y^2 + z^2);
    th = atan2(y, z); 
    rdot = (y*ydot + z*zdot) / r;
    thdot = (y*zdot - z*ydot) / r^2;

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
