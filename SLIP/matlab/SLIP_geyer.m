%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spinrg Loaded Inverted Pendulum (SLIP), Geyer's model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP params
params.m = 22;           % CoM mass (Achilles mass 22 kg)
params.g = 9.81;         % gravity
params.l0 = 0.5;         % spring free length (Achilles leg length 0.7 m)
params.k = 15000;        % spring stiffness (Achilles 82000 N/m)
params.K = [0.05, 0.17]; % Raibert controller gains
params.dt = 0.001;       % time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial conditions
x0_cart = [0.0;   % px
           1.0;   % pz
           1.0;   % vx
           1.0];  % vz

% inital impact angle
alpha = angle_control(x0_cart, params);

% flight phase
tic;
[tspan_f, x_flight] = SLIP_flight(x0_cart, params);
msg = ['Time to compute flight phase dynamics: ', num2str(toc), ' seconds'];

% compute the foot position at impact
px_foot = x_flight(end,1) + params.l0 * cos(alpha);

% convert to polar coordinates
x0_flight = x_flight(end,:);
x0_polar = cart_to_polar(x0_flight, alpha, params);

% forward prop the ground phase
[tspan_g, x_ground] = SLIP_ground(x0_polar, params);

% conert the polar coordinaates to cartesian
for i = 1:length(tspan_g)
    x_cart = polar_to_cart(x_ground(i,:), params);
    x_ground(i,:) = x_cart;
    x_ground(i,1) = x_ground(i,1) + px_foot;
end

% plot
X = [x_flight; x_ground];
plot(X(:,1), X(:,2), 'k', 'LineWidth', 2);
xlabel('x'); ylabel('z');
yline(0)
title('SLIP Geyer Model');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEYER'S MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% closed form solution of the ground dynamics
function [tspan, x_polar] = SLIP_ground(x0_polar, params)

    % unpack params
    m = params.m;
    g = params.g;
    l0 = params.l0;
    k = params.k;
    dt = params.dt;

    % unpack the initial condition, x0
    r = x0_polar(1);
    phi = x0_polar(2);
    r_dot = x0_polar(3);
    phi_dot = x0_polar(4);
    alpha = pi - phi;

    % compute the constant Energy (E) and angular momentum (P)
    P = m * r^2 * phi_dot;
    E = 0.5 * m * r_dot^2 + 0.5 * P^2 / (m * r^2) + 0.5 * k * (l0 - r)^2 + m * g * r;
    epsilon = 2 * E / (m * l0^2);
    omega = P / (m * l0^2);
    omega0 = sqrt(k/m);
    b = sqrt((omega^2 - g/l0)^2 + (omega0^2 + 3*omega^2) * (epsilon - omega^2 - 2*g/l0)) / (omega0^2 + 3*omega^2);

    % compute constant terms
    omega_hat = sqrt(k/m + 3 * phi_dot^2);

    % compute the total stance phase time
    ts = (1  /omega_hat) * (pi + 2 * atan2((g - l0*phi_dot^2), (abs(r_dot)*omega_hat)));
    tspan = linspace(0, ts, ts/dt + 1);

    % get the radial solution, r(t) and phi(t)
    R = zeros(length(tspan),1);
    Phi = zeros(length(tspan),1);
    Rdot = zeros(length(tspan),1);
    Phidot = zeros(length(tspan),1);
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

        % compute the radial rate (Differentiate Eq. 12)
        rdot_t = l0 * b * cos(omega_hat * t) * omega_hat;
        Rdot(i) = -rdot_t;

        % compute the angular rate (Eq. 16)
        rho_t = (R(i) - l0) / l0;
        phi_dot = omega / (1 + rho_t)^2;
        Phidot(i) = phi_dot;
    end

    % pack the polar state
    x_polar = [R, Phi, Rdot, Phidot];
end

% closed form solution of the flight dynamics
function [tspan, x_cart] = SLIP_flight(x0_cart, params)

    % unpack params
    g = params.g;
    l0 = params.l0;
    dt = params.dt;

    % extract inital condition parameters
    x0 = x0_cart(1);
    z0 = x0_cart(2);
    xdot0 = x0_cart(3);
    zdot0 = x0_cart(4);

    % compute the impact angle
    alpha = angle_control(x0_cart, params)

    % compute the flight time
    zf = l0 * sin(alpha);
    tf_sols = roots([-0.5 * g, zdot0, z0 - zf]);
    tf = max(tf_sols);
    tspan = linspace(0, tf, tf/dt + 1);

    % compute the projectile motion
    px = zeros(length(tspan),1);
    pz = zeros(length(tspan),1);
    vx = zeros(length(tspan),1);
    vz = zeros(length(tspan),1);
    for i = 1:length(tspan)
        
        % get the current time
        t = tspan(i);

        % compute the positions
        x_t = x0 + xdot0 * t;
        z_t = z0 + zdot0 * t - 0.5 * g * t^2;
        px(i) = x_t;
        pz(i) = z_t;

        % compute the velocities
        vx_t = xdot0;
        vz_t = zdot0 - g * t;
        vx(i) = vx_t;
        vz(i) = vz_t;
    end

    % pack the state
    x_cart = [px, pz, vx, vz];
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
    alpha = K_p * (px_actual - px_des) + K_v * (vx_actual - vx_des) + pi/2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARTESIAN <---> POLAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cartesian to polar (Eq. 24)
function x_polar = cart_to_polar(x_cart, alpha, params)
    
    % unpack the cartesian state
    vx = x_cart(3);
    vz = x_cart(4);

    % convert to polar
    r = params.l0;
    phi = pi - alpha;
    r_dot = vx * cos(phi) + vz * sin(phi);
    phi_dot = (vz * cos(phi) - vx * sin(phi)) / params.l0;

    % pack the polar state
    x_polar = [r; phi; r_dot; phi_dot];
end

% polar to cartesian (did this myself)
function x_cart = polar_to_cart(x_polar, params)
    
    % unpack the polar state
    r = x_polar(1);
    phi = x_polar(2);
    r_dot = x_polar(3);
    phi_dot = x_polar(4);

    % convert to cartesian
    px = r * cos(phi);
    pz = r * sin(phi);
    vx = r_dot * cos(phi) - r * phi_dot * sin(phi);
    vz = r_dot * sin(phi) + r * phi_dot * cos(phi);

    x_cart = [px, pz, vx, vz];
end
