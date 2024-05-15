%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability of the HLIP Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% HLIP Parameters
T_DSP = 0;
T_SSP = 0.15;
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

% LQR gains
Q = diag([1, 1]);
R = 1;
K_lqr = dlqr(A, B, Q, R);

% closed loop A matrix
A_cl = A - B * K_lqr;

%------------------------------ Define Sets ---------------------

% state constraint set
H_x = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_x = 1 * [1; 1; 1; 1];

% input constraint set
u_lim = 0.1;
H_u = [1; 
       -1];
h_u = u_lim * [1; 1];
U = u_lim * [-1; 1];

% inital set X0, forwards reachability
H_x0 = [1, 0; 
        0, 1; 
       -1, 0;
        0, -1];
h_x0 = 0.1 * [1; 1; 1; 1];

% target set S, backwards reachability
H_s = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_s = 0.1 * [1; 1; 1; 1];

%------------------------------ Backwards Reachability ---------------------
figure(); hold on; grid on; axis equal;

% compute the first precursor set
H  = [H_s * A, H_s * B;
      zeros(2,2), H_u];
h = [h_s; h_u];
Pre_S = con2vert(H, h);           % <-- Why does this work sometimes?
Pre_S = Pre_S(convhull(Pre_S),:);

trisurf(convhull(Pre_S), Pre_S(:,1), Pre_S(:,2), Pre_S(:,3), 'FaceColor', 'b', 'FaceAlpha', 0.5);
legend(["$Pre(S) \times \mathcal U$"], 'interpreter', 'latex');
xlabel("$x_1$", 'interpreter', 'latex');
ylabel("$x_2$", 'interpreter', 'latex');
zlabel("$u$", 'interpreter', 'latex');

% % ok start plotting in 2D
% figure(); hold on; grid on; axis equal;

% % plot the target set
% V_S = con2vert(H_s, h_s);
% V_S = V_S(convhull(V_S),:);
% S = polyshape(V_S);
% plot(S);

% % find the first backward reachable set
% V_Xcal = Pre_S(:,1:2);
% V_Xcal = V_Xcal(convhull(V_Xcal),:)
% Pre_S_Xcal = polyshape(V_Xcal);
% plot(Pre_S_Xcal);

%------------------------------ Forwards Reachability ---------------------
figure(); hold on; grid on; axis equal;

% initial state set
X0 = con2vert(H_x0, h_x0);
X0 = X0(convhull(X0),:);
X0_ = unique(X0, 'rows');
p_X0 = polyshape(X0);
plot(p_X0);

% minkoswki sum (do it this way because U is one-dim)
AX = X0_ * A';
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU); 
plot(p_AX_sum_BU);

% do it again!
AX = AX_sum_BU * A';
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% last time, I swear!
AX = AX_sum_BU * A';
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

legend(["$X_0$", "$Suc(X_0)$", "$Suc(Suc(X_0))$", "$Suc(Suc(Suc(X_0)))$"], 'interpreter', 'latex');

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