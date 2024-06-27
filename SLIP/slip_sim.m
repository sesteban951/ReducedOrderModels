%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple SLIP simualation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% system parameters
params.g = -9.81;    % gravity
params.l = 1.0;     % leg length
params.m = 1.0;     % mass
params.k = 20.0;    % spring constant
params.K = 0.0;     % controller gain

% simulation times
tspan = [0, 5];

% % intial flight condition
x0 = [0;   % y
      3;   % z
      0.1; % ydot
      0];  % zdot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPORT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SLIP ground dynamics
function xdot = slip_ground_dynamics(x_polar, params)

    % unpack the states
    r = x_polar(1);
    th = x_polar(2);
    rdot = x_polar(3);
    thdot = x_polar(4);

    rddot = r * thdot^2 - params.g * cos(th) + (params.k/params.m) * (params.l - r);
    thddot = -2 * rdot * thdot / r - params.g * sin(th) / r;

    % build the state vector
    xdot = [rdot; 
            thdot; 
            rddot; 
            thddot];
end

% get the SLIP ballistic trajecotry
function [tspan, x] = slip_flight_trajectory(x0_cart, alpha, params)
    
    % unpack the intial cartesian state
    z0 = x0_cart(2);
    zdot0 = x0_cart(4);

    % find the leg offset
    pz_touchdown = params.l * cos(alpha);
    
    % find the flight time where the leg touches down
    t_touchdown = roots([-0.5 * params.g, zdot0, z0 - pz_touchdown]);
    t_touchdown = max(t_touchdown);

    % create a tspan upto an including touchdown
    tspan = linspace(0, t_touchdown, 100);

    % compute the ballistic trajectory
    x = zeros(4, length(tspan));
    for i = 1 : length(tspan)
        x_t = slip_flight_sol(tspan(i), x0_cart, params);
        x(:, i) = x_t';
    end
end

% SLIP flight dynamics
function x_t = slip_flight_sol(t, x0_cart, params)
    
    % unpack the states
    y0 = x0_cart(1);
    z0 = x0_cart(2);
    ydot0 = x0_cart(3);
    zdot0 = x0_cart(4);

    % linear ballistic dynamics
    yf = y0 + ydot0 * t;
    zf = z0 + zdot0 * t - 0.5 * params.g * t^2;
    ydotf = ydot0;
    zdotf = zdot0 - params.g * t;

    % solution to linear dynamics
    x_t = [yf; zf; ydotf; zdotf];
end

% raibert foot placement controller
function alpha = raibert_controller(x_cart, params)
    
    % unpack the states
    ydot = x_cart(3);

    % raibert controller
    alpha = params.K * (ydot - 0.0);
end

% convert all cartesian states to polar states
function X_polar = convert_cart_to_polar_states(x_cart)
    
    % get how many states we have
    [r, c] = size(x_cart);

    % convert every state into cartesian
    X_polar = zeros(r, c);
    for i = 1 : r
        x_polar = cart_to_polar(x_cart(i, :));
        X_polar(i, :) = x_polar';
    end
end

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

% convert all polar states to cartesian states
function X_cart = convert_polar_to_cart_states(x_polar)
    
    % get how many states we have
    [r, c] = size(x_polar);

    % convert every state into cartesian
    X_cart = zeros(r, c);
    for i = 1 : r
        x_cart = polar_to_cart(x_polar(i, :));
        X_cart(i, :) = x_cart';
    end
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