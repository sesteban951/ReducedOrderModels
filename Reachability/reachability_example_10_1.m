%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.1, discrete-time linear system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k
A = [0.5, 0;
     1 -0.5];

% Xcal constraint set. Convert to H-polytope, S = {x | Hx <= h}
H = [1, 0; 
     0, 1; 
    -1, 0;
     0, -1];
h = [10; 10; 10; 10];

% ---------------------- Bacwards Reachability -----------------------------
% Pre(S) is the set of states which can be driven into the target set S in
% one time step.

% target set, S
S = get_polyshape(H, h);

% one step backwards reachable set, Pre(Xcal)
C = H * A;
d = h;
Pre_S = get_polyshape(C, d);
figure; hold on;
plot(S);
plot(Pre_S);
legend(["$S$", "$Pre(S)$"], 'interpreter', 'latex');
grid on; axis equal;

% one step constrained backwards reachable set, Pre(Xcal) intersect Xcal
C = [H * A; H];
d = [h; h];
Pre_S_S = get_polyshape(C, d);
figure; hold on;
plot(S);
plot(Pre_S_S);
legend(["$S$", "$Pre(S) \cap S$"], 'interpreter', 'latex');
grid on; axis equal;

% ---------------------- Forward Reachability -----------------------------
% Every state in the set X0 can be driven into the resulting set, Suc(X0) 
% in one time step.

% initial state set
X0 = con2vert(H, h);
k = convhull(X0);
V = X0(k,:);
% [V, v] = vert2con(V); % just to check that the same set is obtained, looks good

% one step forward reachable set, Suc(X0)
AV = V * A';
k = convhull(AV);
V = AV(k,:);
[C, d] = vert2con(V);
Suc_X0 = get_polyshape(C, d);
figure; hold on;
plot(get_polyshape(H, h));
plot(Suc_X0);
legend(["$X_0$", "$Suc(X_0)$"], 'interpreter', 'latex');
grid on; axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = get_polyshape(A, b)

     % convert to V-polytope, this is an addon function
     V = con2vert(A, b);

     % get convex hull so no repeated vertices, Recall that
     % image of a convex set under an affine transformation is convex
     k = convhull(V);
     V = V(k,:);

     % create polyshape
     p = polyshape(V, 'Simplify', true);
end
