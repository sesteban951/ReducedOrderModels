%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.4, "Acutauted N-step Forwards Reachable Set"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k
A = [1.5, 0;
     1 -1.5];
B = [1; 0];

% intial state set, X0 = {x0 | Hx0 <= h0}
H_x0 = [1, 0; 
        0, 1; 
       -1, 0;
        0, -1];
h_x0 = [1; 1; 1; 1];

% state constraint set, Xcal = {x | Hx <= h}
H_x = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_x = [10; 10; 10; 10];

% input constraint set, U = {u | Hu <= h}
H_u = [1; 
      -1];
h_u = [5; 5];
U = [-5; 5];

% ---------------------- Forwards Reachability -----------------------------
% I COULDN'T GET THE RESULTS IN THE BOOK.

figure(); grid on; hold on;

% plot the orignal intial condition set
X0 = con2vert(H_x0, h_x0);
X0 = X0(convhull(X0),:);
p_X0 = polyshape(X0);
plot(p_X0);

% plot the one step reachable set
AX = X0 * A'
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% do another one
AX = AX_sum_BU * A';
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];     
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% do another one
AX = AX_sum_BU * A';
BU = U * B';
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];     
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);
