%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability of the HLIP Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% HLIP Parameters
T_DSP = 0;
T_SSP = 0.25;
g = 9.81;
z0 = 0.5;
lambda = sqrt(g/z0);
A_SSP = [0, 1;
         lambda^2, 0];

% step-to-step discrete dynamics
A = exp(A_SSP * T_SSP) * [1, T_DSP; 0, 1];
B = exp(A_SSP * T_SSP) * [-1; 0];

% disturbance limits
w1_lims = [-1, 1];
w2_lims = [-1, 1];

% Deadbeat gains
coth_ = coth(T_SSP * lambda);
Kp_db = 1;
Kd_db = T_DSP + (1/lambda) * coth_;
K_db = [Kp_db, Kd_db];

% LQR gains
Q = diag([1, 1]);
R = 1;
K_lqr = dlqr(A, B, Q, R);

% closed loop A matrix
A_cl_db = A + B * K_db;
eig(A_cl_db);
A_cl_lqr = A - B * K_lqr;
eig(A_cl_lqr);

% initial state
x0 = [1; 5];

% forward propagation
N = 2;
x = forward_propagate(A_cl_lqr, x0, N);

% plot the results
figure(1); hold on; grid on;
xline(0); yline(0);
plot(x(:, 1), x(:, 2), 'ko', 'LineWidth', 2);
plot(x0(1), x0(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x(end, 1), x(end, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('$p$', 'Interpreter', 'latex');
ylabel('$v$', 'Interpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxilllary tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% forward propagation of the HLIP model for N steps
function X = forward_propagate(A_cl, x0, N)
    
    % initialize state vector history
    X = zeros(N+1, 2);

    % initial state
    X(1, :) = x0';

    % forward propagation
    x_k = x0;
    for i = 1:N
        x_k = A_cl * x_k;
        X(i+1, :) = x_k';
    end
    
end

% sample random distruabnce
function w = rand_w(w1_lims, w2_lims)


    % unform sample of the disturbance w
    w = [x1_lims(1) + (x1_lims(2) - x1_lims(1)) * rand();
         x2_lims(1) + (x2_lims(2) - x2_lims(1)) * rand()];
end

