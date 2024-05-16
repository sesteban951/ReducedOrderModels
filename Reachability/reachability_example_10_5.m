%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.5, "Autonomous Control Invariant Set"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k
A = [0.5, 0;
     1 -0.5];

% state constraint set
H_x = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_x = 10 * [1; 1; 1; 1];

% plot the state constraint set
figure(); grid on; axis equal; hold on;
X = con2vert(H_x, h_x);
X = X(convhull(X),:);
p_X = polyshape(X);
plot(p_X);

% plot the one step backwards reachable set
C = [H_x * A; H_x];
d = [h_x; h_x];
Pre_X = con2vert(C, d);
Pre_X = Pre_X(convhull(Pre_X),:);
p_Pre_X = polyshape(Pre_X);
plot(p_Pre_X);

% plot another one (TURNS OUT YOU CONVERGE IN ONE ITERATION)
C = [C * A; H_x];
d = [d; h_x];
Pre_X = con2vert(C, d);
Pre_X = Pre_X(convhull(Pre_X),:);
p_Pre_X = polyshape(Pre_X);
plot(p_Pre_X);

legend(["$X$", "$Pre(X)$", "$Pre(Pre(X))$"], 'interpreter', 'latex');

% print the O_inf values (looks right)
O = con2vert(C, d);
O = O(convhull(O),:);
[O_inf, o_inf] = vert2con(O)

