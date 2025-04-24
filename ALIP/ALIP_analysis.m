%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALIP simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system parameters
params.m = 5.0;   % mass [kg]
params.g = 9.81;  % gravity [m/s^2]
params.z0 = 1.0;  % center of mass height [m]

% simulation params
params.T = 0.5;  % time to complete a step [sec]

% foot placement gains
params.kp_fp = 1.0; % foot placement position gain
params.kd_fp = 0.5; % foot placement velocity gain
params.l_max = 1.0; % maximum foot placement [m]

% ankle torque gains
params.kp_ankle = 10.0; % ankle torque position gain
params.kd_ankle = 1.0;  % ankle torque velocity gain
params.tau_max = 3.0; % maximum ankle torque [Nm]

% input options
params.use_capture_point = 1; % use capture point for foot placement (0 or 1)
params.use_ankle_torque = 0;  % use ankle torque for control (0 or 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial state
t0 = 0;
x0 = [0.25;       % p, position of the COM relative to stance
      0.0];      % L, angular momentum about stance
p_stance = 0.0;  % initial position of the stance foot

% simulation time
tspan = 0: 0.01 : params.T;
n_steps = 10;

% forward simulation
t_sol = zeros(length(tspan), n_steps);
x_sol = zeros(length(tspan), 2, n_steps);
p_stance_sol = zeros(length(tspan), n_steps);
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

    % compute the new stance foot position
    p_stance_sol(:, i) = p_stance;
    p_stance = p_stance + l;

    % apply the reset map
    x0 = reset_map(x_minus, l);
end

% stack all the COM states and time
p_com_W = [];
p_stance_W = [];
t = [];
for i = 1:n_steps
    p_com = x_sol(:, 1, i) + p_stance_sol(:, i);
    p_stance = [p_stance; p_stance_sol(:, i)];
    t_series = t_sol(:, i);

    p_com_W = [p_com_W; p_com];
    p_stance_W = [p_stance_W; p_stance_sol(:, i)];
    t = [t; t_series];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); 
set(gcf, 'Position', get(0, 'ScreenSize')); % Set figure to full screen 


% plot the position of the COM
subplot(2, 2, 1);
hold on;
for i=1:n_steps
    plot(t_sol(:, i), x_sol(:, 1, i), 'b', 'LineWidth', 2);
end
yline(0);
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
yline(0);
xlabel('Time (s)');
ylabel('L');
title('Angular momentum about stance');
grid on;

% plot the state trajectory
subplot(2, 2, 2);
hold on; 
grid on;
for i=1:n_steps
    plot(x_sol(:, 1, i), x_sol(:, 2, i), 'k', 'LineWidth', 2);
end
xline(0); yline(0);
xlabel('p [m]');
ylabel('L [kg m/s]');
title('State trajectory');

% animate the ALIP
subplot(2, 2, 4);
hold on;
xline(0); yline(0);
grid on; axis equal;

% find the plot limits
x_min = min([p_com_W; p_stance_W]) - 0.25;
x_max = max([p_com_W; p_stance_W]) + 0.25;
y_min = -0.1;
y_max = params.z0 + 0.25;
xlim([x_min, x_max]);
ylim([y_min, y_max]);

tic; 
ind = 1;
while ind < length(t)

    % current COM position
    p_com = p_com_W(ind, :);

    % current foot position
    p_foot = p_stance_W(ind, :);

    % plot the ALIP
    pole = plot([p_foot, p_com], [0, params.z0], 'k', 'LineWidth', 2);
    mass = plot(p_com, params.z0, 'ko', 'MarkerSize', 30, 'MarkerFaceColor', "#D95319", 'LineWidth', 2);
    foot = plot(p_foot, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

    drawnow;

    % set the title
    title(['t = ', num2str(t(ind), '%.2f'), ' s']);

    % update the index
    ind = ind + 1;
    while toc < t(ind)
        % wait until time to update
    end

    % remove the previous plot
    if ind ~= length(t)
        delete(pole);
        delete(mass);
        delete(foot);
    end
end

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

    % saturate the foot placement
    if l > params.l_max
        l = params.l_max;
    elseif l < -params.l_max
        l = -params.l_max;
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
