%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CART POLE PLotting the Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% import the saved data
state_data = importdata('data/state_data.csv');
control_data = importdata('data/input_data.csv');
time_data = importdata('data/time_data.csv');

% extract the time data
time = time_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_states = 0;
animation = 1;

realtime_rate = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some colors
gray = [0.5, 0.5, 0.5];
orange = [1, 0.5, 0];

% length of the pole
l = 0.5; 

% select time range
% t0 = time_data(1);
% tf = time_data(end);
t0 = 0.0;
tf = 7.8;

% find the indices for the selected time range
idx = find(time >= t0 & time <= tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the time
time = time(idx);

% extract the states
pos = state_data(:,1);   % cart pole position
pos = pos(idx);
theta = state_data(:,2); % pole angle
theta = theta(idx);
vel = state_data(:,3);   % cart velocity
vel = vel(idx);
omega = state_data(:,4); % pole angular velocity
omega = omega(idx);

% extract the control inputs
force = control_data(idx(1:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_states == 1
    % plot all the data
    figure('Name','Cart Pole Trajectory Optimization Results');

    % Cart Position
    subplot(2,3,1);
    plot(time, pos, 'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Cart Position [m]');
    grid on;

    % Pole Angle
    subplot(2,3,2);
    plot(time, theta, 'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Pole Angle [rad]');
    grid on;

    % Cart Velocity
    subplot(2,3,4);
    plot(time, vel, 'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Cart Velocity [m/s]');
    grid on;

    % Pole Angular Velocity
    subplot(2,3,5);
    plot(time, omega, 'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Pole Angular Velocity [rad/s]');
    grid on;

    % Control Input
    subplot(2,3,[3, 6]);
    plot(time, force, 'r','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Force [N]');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if animation == 1
    figure('Name','Cart Pole Animation');
    hold on; grid on; axis equal;
    xline(0); yline(0);
    xlabel('x [m]'); ylabel('y [m]');
    x_lims = [min(pos) - 0.25, max(pos) + 0.25];
    y_lims = [-l-0.1, l + 0.2];
    xlim(x_lims); ylim(y_lims);

    % scale the time date to speed up or slow down the animation
    time = time / realtime_rate;

    tic;
    t_now = time(1);
    ind = 1;
    while t_now < time(end)

        % get the current position
        p = pos(ind);
        th = theta(ind);

        % cart and pole positions
        cart_pos = [p; 0];
        pole_pos = [p + l*sin(th); -l*cos(th)];

        % draw the cart
        cart = rectangle('Position', [cart_pos(1)-0.1, cart_pos(2)-0.05, 0.2, 0.1], 'Curvature', 0.1, 'FaceColor', gray);
        pole_base = plot(cart_pos(1), cart_pos(2), 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k');

        % draw pole
        pole = plot([cart_pos(1), pole_pos(1)], [cart_pos(2), pole_pos(2)], 'k', 'LineWidth', 2);
        ball = plot(pole_pos(1), pole_pos(2), 'ko', 'MarkerSize', 16, 'MarkerFaceColor', orange);

        % show the current time
        msg = sprintf('Time: %.2f s', time(ind) * realtime_rate);
        title(msg);

        % update the plot
        drawnow;

        % wait unitl the next time step
        while toc < time(ind+1)
            % wait
        end

        % increment the index
        if ind+1 == length(time)
            break;
        else
            ind = ind + 1;

            delete(cart)
            delete(pole_base)
            delete(pole)
            delete(ball)
        end
    end

end
