%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.9, "Actuated system with uncertainty and inputs"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k + w_k
A = [1.5, 0;
     1 -1.5];
B = [1; 0];

% input constraint set, Ucal
H_u = [1; -1];
h_u = 5.0 * [1; 1];
U = [-5; 5];

% state constraint set, Xcal
H_x = [1, 0;
       0, 1;
      -1, 0;
       0, -1];
h_x = 10 * [1; 1; 1; 1];

% distrubance set, Wcal
H_w = [1, 0'
       0, 1;
      -1, 0;
       0, -1];
h_w = 1.0 * [1; 1; 1; 1];

% -------------------------- Forward Reachability ---------------------------
figure(); hold on; grid on; axis equal;
X = con2vert(H_x, h_x);
X = X(convhull(X),:);
p_X = polyshape(X);
plot(p_X);

% Compute one step Pre(X, U, W), Eq. (10.49)
H = [H_x * A, H_x * B;
     zeros(size(H_u, 1), size(H_x, 2)), H_u];

% solve the linear program (Eq. 10.50)
% h_i_tilde = min (h_i - H_i * w) subject to  H_w * w <= h_w
h_tilde = zeros(length(h_x), 1);
for i = 1:length(h_x)

     H_i = H_x(i,:);
     h_i = h_x(i);

     [w, cost] = linprog(-H_i, H_w, h_w);
     h_tilde(i) = h_i - H_i * w;
end
h = [h_tilde; h_u];

% plot Pre(X, U, W) 
V = con2vert(H, h);
V_x = V(:,1:2);
V_x = V_x(convhull(V_x),:);
p_V_x = polyshape(V_x);
plot(p_V_x);

% plot Pre(X, U, W) intersect X
H_ = [H;
      H_x, zeros(size(H_x, 1), size(H_u, 2))];
h_ = [h; h_x];

V = con2vert(H_, h_);
V_x = V(:,1:2);
V_x = V_x(convhull(V_x),:);
p_V_x = polyshape(V_x);
plot(p_V_x);

legend(["$\mathcal{X}$", "$Suc(\mathcal{X}, \mathcal{U}, \mathcal{W})$", "$Suc(\mathcal{X}, \mathcal{U}, \mathcal{W}) \cap \mathcal{X}$"], 'interpreter', 'latex');


% compute  Suc(X, U, W)
AX = X * A';
AX = unique(AX, 'rows');
W = con2vert(H_w, h_w);

AX_sum_W = [];
for i = 1:size(W, 1)
    w = W(i,:);
    
    AX_sum_W = [AX_sum_W; AX + ones(size(AX,1),length(w)) * w'];
end

BU = U * B';
AX_sum_W_sum_BU = [];
for i = 1:size(BU, 1)

     bu = BU(i,:);
     AX_sum_W_sum_BU = [AX_sum_W_sum_BU; AX_sum_W + ones(size(AX_sum_W,1),length(bu)) * bu'];
end  

AX_sum_W_sum_BU = AX_sum_W_sum_BU(convhull(AX_sum_W_sum_BU),:);
p_AX_sum_W_sum_BU = polyshape(AX_sum_W_sum_BU);
plot(p_AX_sum_W_sum_BU);


