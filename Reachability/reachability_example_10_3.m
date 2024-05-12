%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.3, discrete-time linear system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k
A = [1.5, 0;
     1 -1.5];
B = [1; 0];

% target set, S = {s | Hs <= h}
H_s = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_s = [1; 1; 1; 1];

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

% ---------------------- Backwards Reachability -----------------------------

% compute Pre(S) x U
H = [H_s * A,    H_s * B;
     zeros(2,2), H_u];
h = [h_s; h_u];
Pre_S_V = con2vert(H, h);

% plot 3D polytope in Scal x Ucal
figure(); grid on; axis equal;
trisurf(convhull(Pre_S_V), Pre_S_V(:,1), Pre_S_V(:,2), Pre_S_V(:,3), 'FaceColor', 'b', 'FaceAlpha', 0.5);
legend(["$Pre(S) \times \mathcal U$"], 'interpreter', 'latex');
xlabel("$s_1$", 'interpreter', 'latex');
ylabel("$s_2$", 'interpreter', 'latex');
zlabel("$u$", 'interpreter', 'latex');

% plot Scal;
V = con2vert(H_s, h_s);
V = V(convhull(V),:);
figure(); grid on; axis equal; hold on;
plot(V(:,1), V(:,2), '-k', 'LineWidth', 2); 

% plot 2D projection onto Xcal.
Pre_S_V_Xcal = Pre_S_V(:,1:2)
Pre_S_V_Xcal = Pre_S_V_Xcal(convhull(Pre_S_V_Xcal),:)
plot(Pre_S_V_Xcal(:,1), Pre_S_V_Xcal(:,2), 'r', 'LineWidth', 2);

[H_, h_] = vert2con(Pre_S_V_Xcal);
C = [H_ * A, H_ * B; 
     zeros(2,2), H_u];
d = [h_; h_u];
V = con2vert(C, d);
V = V(:,1:2);
V = V(convhull(V),:);
plot(V(:,1), V(:,2), 'g', 'LineWidth', 2);

Pre_S_V_Xcal = V;

[H_, h_] = vert2con(Pre_S_V_Xcal);
C = [H_ * A, H_ * B; 
     zeros(2,2), H_u];
d = [h_; h_u];
V = con2vert(C, d);
V = V(:,1:2);
V = V(convhull(V),:);
plot(V(:,1), V(:,2), 'cyan', 'LineWidth', 2);

Pre_S_V_Xcal = V;

[H_, h_] = vert2con(Pre_S_V_Xcal);
C = [H_ * A, H_ * B; 
     zeros(2,2), H_u];
d = [h_; h_u];
V = con2vert(C, d);
V = V(:,1:2);
V = V(convhull(V),:);
plot(V(:,1), V(:,2), 'k', 'LineWidth', 2);

legend(["$\mathcal{S}$", "$\mathcal{B}_1(S)$", "$\mathcal{B}_2(S)$", "$\mathcal{B}_3(S)$"], 'interpreter', 'latex');

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