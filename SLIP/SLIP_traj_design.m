%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP Traj Design (based on paper Wensing and Orin 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;



v_min = 0.0;
v_max = 5.0;

v = v_min: 0.05 :v_max;

% plot the cadence
c = cadence(v);

% plot the desired stance time
t_s = stance_time(v);
plot(v, t_s, 'r', 'LineWidth', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cadence parameters
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