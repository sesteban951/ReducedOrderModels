%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.2, discrete-time linear system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k
A = [1.5, 0;
     1 -1.5];
B = [1; 0];

% target set, S = {x | Hx <= h}
H_x = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_x = [10; 10; 10; 10];

% target set, U = {u | Hu <= h}
H_u = [1; 
      -1];
h_u = [5; 5];

% ---------------------- Backwards Reachability -----------------------------
% Pre(S) is the set of states which can be driven into the target set S in
% one time step. Defines polytope in Xcal x Ucal
H = [H_x * A,    H_x * B;
     zeros(2,2), H_u];
h = [h_x; h_u];
V = con2vert(H, h);

% plot 3D polytope in Xcal x Ucal
figure(); grid on; axis equal;
trisurf(convhull(V), V(:,1), V(:,2), V(:,3), 'FaceColor', 'b', 'FaceAlpha', 0.5);
legend(["$Pre(S) \times \mathcal U$"], 'interpreter', 'latex');
xlabel("$x_1$", 'interpreter', 'latex');
ylabel("$x_2$", 'interpreter', 'latex');
zlabel("$u$", 'interpreter', 'latex');

% plot 2D projection onto Xcal. Recall that projection of a convex set
% to one of its coordinates is convex (Boyd ch. 2.3.2)
V_x = con2vert(H_x, h_x);
V_x = V_x(convhull(V_x),:);
S = polyshape(V_x);

V_Xcal = V(:,1:2);
k = convhull(V_Xcal);
V_Xcal = V_Xcal(k,:);
Pre_S = polyshape(V_Xcal);

[C_Xcal, d_Xcal] = vert2con(V_Xcal);
C = [C_Xcal; H_x];
d = [d_Xcal; h_x];
Pre_S_S = get_polyshape(C, d);

% plot 2D polytopes in Xcal
figure(); grid on; axis equal; hold on;
plot(S);
plot(Pre_S);
plot(Pre_S_S);
legend(["$S$", "$Pre(S)$", "$Pre(S) \cap S$"], 'interpreter', 'latex');

% ---------------------- Forwards Reachability -----------------------------

X0 = con2vert(H_x, h_x);
k = convhull(X0);
X0 = X0(k,:);
X0 = remove_redundant(X0);

% (A o X) = conv(AV)
AX = X0 * A'

% (B o U) = conv(BU)
U = [-5; 5];
BU = U * B';

% minkowsli sum of AX and BU
AX_sum_BU = [AX + ones(size(AX,1),1) * BU(1,:);
             AX + ones(size(AX,1),1) * BU(2,:)];
k = convhull(AX_sum_BU);
AX_sum_BU = AX_sum_BU(k,:);
AX_sum_BU = polyshape(AX_sum_BU);

% plot 2D polytopes in Xcal
figure(); grid on; axis equal; hold on;
plot(get_polyshape(H_x, h_x));
plot(AX_sum_BU);
legend(["$\mathcal{X}_0$", "$Succ(\mathcal{X}_0)$"], 'interpreter', 'latex');


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

function v = remove_redundant(V)
    
    % get convex hull, may have repeated vertices
    k = convhull(V);
    V = V(k,:);

    % get unique vertices
    v = unique(V, 'rows');

end