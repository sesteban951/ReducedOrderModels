%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.8, "Autonomous system with uncertainty"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k + w_k
A = [0.5, 0;
     1 -0.5];

% state constraint set, Xcal
H_x = [1, 0;
       0, 1;
      -1, 0;
       0, -1];
h_x = 10 * [1; 1; 1; 1];

% distrubance set, Wcal
H_w = [1, 0;
       0, 1;
      -1, 0;
       0, -1];
h_w = 1.0 * [1; 1; 1; 1];

% ------------------------ Backwards Reachability -----------------------------

figure(); hold on; grid on; axis equal;
X = con2vert(H_x, h_x);
X = X(convhull(X),:);
p_X = polyshape(X);
plot(p_X);

% Eq. (10.38) and (10.39)
H = H_x * A;

% solve the linear program (Eq. 10.40)
% h_i_tilde = min (h_i - H_i * w) subject to  H_w * w <= h_w
h_tilde = zeros(length(h_x), 1);
for i = 1:length(h_x)

     H_i = H_x(i,:);
     h_i = h_x(i);
     
     [w, cost] = linprog(-H_i, H_w, h_w);
     h_tilde(i) = h_i - H_i * w;
end

% compute Pre(X, W) = {x | H_x * A <= h_x - H_x * w}
H = [H_x * A; H_x];
h = [h_tilde; h_x];

V = con2vert(H, h);
V = V(convhull(V),:);
p_V = polyshape(V);
plot(p_V);

legend(["$\mathcal{X}$", "$Pre(\mathcal{X}, \mathcal{W})$"], 'interpreter', 'latex');

[C, d] = vert2con(V)
