%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALIP simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system parameters
params.m = 1.0;   % mass
params.g = 9.81;  % gravity
params.z0 = 1.0;  % center of mass height

% simulation params
params.T = 0.4;  % time to complete a step

% foot placement gains
params.kp_fp = 1.0; % foot placement position gain
params.kd_fp = 0.5; % foot placement velocity gain

% ankle torque gains
params.kp_ankle = 10.0; % ankle torque position gain
params.kd_ankle = 1.0;  % ankle torque velocity gain
params.tau_max = 10.0; % maximum ankle torque

% input options
params.use_capture_point = 0; % use capture point for foot placement (0 or 1)
params.use_ankle_torque = 0;  % use ankle torque for control (0 or 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial state
t0 = 0;
x0 = [0.5;       % p, position of the COM relative to stance
      0.4];      % L, angular momentum about stance
p_stance = 0.0;  % initial position of the stance foot

% simulation time
tspan = 0: 0.01 : params.T;
n_steps = 10;

% forward simulation
t_sol = zeros(length(tspan), n_steps);
x_sol = zeros(length(tspan), 2, n_steps);
for i = 1:n_steps
    
    % flow the dynamics forward over the step period
    [~, x] = ode45(@(t,x) alip_continuous_dynamics(t, x, params), tspan, x0);

    % save the time ans state trajectory
    t_sol(:, i) = t0 + tspan';
    x_sol(:, :, i) = x;

    % update the time intial condition
    t0 = t_sol(end, i);
    
    % compute the foot placement
    x_minus = x(end, :)';
    l = alip_foot_placement(x_minus, params);

    % apply the reset map
    x0 = reset_map(x_minus, l);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); 

% plot the position of the COM
subplot(2, 2, 1);
hold on;
for i=1:n_steps
    plot(t_sol(:, i), x_sol(:, 1, i), 'b', 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('p');
title('Position of the COM');
grid on;

% plot the angular momentum
subplot(2, 2, 3);
hold on;
for i=1:n_steps
    plot(t_sol(:, i), x_sol(:, 2, i), 'r', 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('L');
title('Angular momentum about stance');
grid on;

% plot the state trajectory
subplot(2, 2, [2, 4]);
hold on; 
axis equal; grid on;
for i=1:n_steps
    plot(x_sol(:, 1, i), x_sol(:, 2, i), 'k', 'LineWidth', 2);
end
xline(0); yline(0);
xlabel('p [m]');
ylabel('L [kg m/s]');
title('State trajectory');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute a solution to the ALIP dynamics
function xdot = alip_continuous_dynamics(t, x, params)

    % ALIP continuous dynamics matrix
    m = params.m;
    g = params.g;
    z0 = params.z0;

    % dynamics matrices
    A = [0,   1/(m*z0);
         m*g, 0];
    B = [0; 1];

    % compute the continuous input
    tau = alip_continuous_control(t, x, params);

    % compute the dynamics
    xdot = A * x + B * tau;
end

% compute the control action
function tau = alip_continuous_control(t, x, params)
    
    % current state
    p = x(1); % position of the COM relative to stance
    L = x(2); % angular momentum about stance
    v = L/(params.m *params.z0); % velocity of the COM relative to stance

    % compute the ankle torque
    if params.use_ankle_torque
        tau = params.kp_ankle * (0 - p) + params.kd_ankle * (0 - v);

        % limit the ankle torque
        if tau > params.tau_max
            tau = params.tau_max;
        elseif tau < -params.tau_max
            tau = -params.tau_max;
        end

    % no ankle torque
    else
        tau = 0;
    end
end

% compute the required foot placement
function l = alip_foot_placement(x_minus, params)
    
    % preimpact state
    p_minus = x_minus(1); % position of the COM relative to stance
    L_minus = x_minus(2); % angular momentum about stance
    v_minus = L_minus/(params.m * params.z0); % velocity of the COM relative to stance

    % use capture point
    if params.use_capture_point
        % compute the capture point
        lambda = sqrt(params.g / params.z0);
        l = p_minus + v_minus / lambda;
    % use simple PD
    else
        l = params.kp_fp * (p_minus - 0.0) + params.kd_fp * (v_minus - 0.0);
    end
end

% apply the reset map to the state
function x_plus = reset_map(x_minus, l)
    
    % current state
    p_minus = x_minus(1); % position of the COM relative to stance
    L_minus = x_minus(2); % angular momentum about stance
    
    % compute the new state after the reset map
    p_plus = p_minus - l;
    L_plus = L_minus;
    
    % return the new state
    x_plus = [p_plus; L_plus];
end
