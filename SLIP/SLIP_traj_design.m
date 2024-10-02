%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP Traj Design (based on paper Wensing and Orin 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_running_params = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the SLIP design parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RUNNING PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_running_params == 1

    % choose range of desired velocities
    v_min = 0.0;
    v_max = 3.0;
    v = v_min: 0.05 :v_max;

    % plot the cadence
    subplot(2,1,1);
    c = cadence(v);
    plot(v, c, 'b', 'LineWidth', 2);
    xlabel('Desired Velocity [m/s]');
    ylabel('Cadence [steps/min]');
    title('Cadence vs. Desired Velocity');

    % plot the desired stance time
    subplot(2,1,2);
    t_s = stance_time(v);
    plot(v, t_s, 'r', 'LineWidth', 2);
    xlabel('Desired Velocity [m/s]');
    ylabel('Stance Time [s]');
    title('Stance Time vs. Desired Velocity');

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
% Apex-to-Apex Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




