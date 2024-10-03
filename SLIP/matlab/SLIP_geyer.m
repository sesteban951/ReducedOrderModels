%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spinrg Loaded Inverted Pendulum (SLIP), Geyer's model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; 

% SLIP params
params.m = 22;           % CoM mass (Achilles mass 22 kg)
params.g = 9.81;         % gravity
params.l0 = 0.5;         % spring free length (Achilles leg length 0.7 m)
params.k = 15000;        % spring stiffness (Achilles 82000 N/m)
alpha_max_deg = 60;      % max foot angle from verticle [deg]
params.alpha_max = alpha_max_deg * (pi/180);  % max foot angle [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial conditions
x0 = [0;     % px
      1.0;   % pz
      0.0;   % vx
      0.0];  % vz

tic;
[tspan, X, Z] = SLIP_flight(x0, params);
toc;
msg = "Time to compute SLIP flight: " + string(toc) + " seconds";

figure;
plot(X,Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEYER'S MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% closed form solution of the ground dynamics
function [tspan, R, Phi] = SLIP_ground(t_abs, x0_polar, params)

    % unpack params
    m = params.m;
    g = params.g;
    l0 = params.l0;
    k = params.k;

    % unpack the initial condition, x0
    phi = x0_polar(2);
    r_dot = x0_polar(3);
    phi_dot = x0_polar(4);
    alpha = pi - phi;

    % compute constant terms
    omega_hat = sqrt(k/m + 3 * phi_dot^2);

    % compute the total stance phase time
    ts = (1  /omega_hat) * (pi + 2 * atan2((g - l0*phi_dot^2), (abs(r_dot)*omega_hat)));
    tspan = 0:0.005:ts;
    tspan(end) = ts;

    % get the radial solution, r(t)
    R = zeros(length(tspan),1);
    Phi = zeros(length(tspan),1);
    for i = 1:length(tspan)
        
        % current time
        t = tspan(i);

        % compute the radius (Eq. 21)
        r_t = l0 - (abs(r_dot) / omega_hat) * sin(omega_hat * t) + ((phi_dot^2 * l0 - g) / omega_hat^2) * (1 - cos(omega_hat*t));
        R(i) = r_t;

        % compute the angle (Eq. 22)
        term1 = 1 - 2 * (phi_dot^2 - g/l0) / omega_hat^2;
        term2 = (phi_dot^2 - g/l0) * sin(omega_hat*t) / omega_hat^2 + abs(r_dot) * (1 - cos(omega_hat*t)) / (omega_hat * l0);
        phi_t = pi - alpha + term1 * phi_dot * t + (2*phi_dot / omega_hat) * term2;
        Phi(i) = phi_t;
    end
end

% closed form solution of the flight dynamics
function [tspan, X, Z] = SLIP_flight(x0_cart, params)

    % unpack params
    g = params.g;
    l0 = params.l0;

    % extract inital condition parameters
    x0 = x0_cart(1);
    z0 = x0_cart(2);
    xdot0 = x0_cart(3);
    zdot0 = x0_cart(4);

    % compute the impact angle
    alpha = angle_control(x0_cart, params);

    % compute the flight time
    zf = l0 * sin(alpha);
    tf_sols = roots([-0.5 * g, zdot0, z0 - zf]);
    tf = max(tf_sols);
    tspan = 0: 0.005 : tf;
    tspan(end) = tf;

    % compute the projectile motion
    X = zeros(length(tspan),1);
    Z = zeros(length(tspan),1);
    for i = 1:length(tspan)
        
        % get the current time
        t = tspan(i);

        % compute the positions
        x_t = x0 + xdot0 * t;
        z_t = z0 + zdot0 * t - 0.5 * g * t^2;

        X(i) = x_t;
        Z(i) = z_t;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control the impact angle
function alpha = angle_control(x_cart, params)
    
    % unpack the desired state
    px_des = 1.0;
    vx_des = 0.0;

    % simple Raibert controller
    K = params.K;
    K_p = K(1);
    K_v = K(2);

    % compute the desired angle
    px_actual = x_cart(1);
    vx_actual = x_cart(3);
    alpha = K_p * (px_actual - px_des) + K_v * (vx_actual - vx_des);

    % clip the angle to a range
    alpha_low = -params.alpha_max;
    alpha_high = params.alpha_max;
    alpha = max(alpha_low, min(alpha_high, alpha));
end
