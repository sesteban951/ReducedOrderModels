%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reachability example 10.6, "Actuated Control Invariant Set"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% System dynamics, x_{k+1} = A x_k + B u_k
A = [1.5, 0;
     1 -1.5];
B = [1; 0];

% state constraint set
H_x = [1, 0; 
       0, 1; 
      -1, 0;
       0, -1];
h_x = 10 * [1; 1; 1; 1];

H_u = [1; 
      -1];
h_u = [5; 5];
U = [-5; 5];

% plot the intial condition set
figure(); grid on; axis equal; hold on;
X0 = con2vert(H_x, h_x);
X0 = X0(convhull(X0),:);
p_X0 = polyshape(X0);
plot(p_X0);

% plot the one step backwards reachable set
C = [H_x * A, H_x * B;
     zeros(2,2), H_u];
d = [h_x; h_u];
C_x = [H_x, zeros(size(H_x,1),1)];
d_x = h_x;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

% ---------------------- Iterate like a caveman -----------------------------

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);


% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% plot another one
V = con2vert(C, d);
V_x = V(:,1:2);
[H_x_, h_x_] = vert2con(V_x);
C = [H_x_ * A, H_x_ * B;
     zeros(2,2), H_u];
d = [h_x_; h_u];
C_x = [H_x_, zeros(size(H_x_,1),1)];
d_x = h_x_;
C = [C; C_x];
d = [d; d_x];

Pre_X0 = con2vert(C, d);
Pre_X0 = Pre_X0(:,1:2);
Pre_X0 = Pre_X0(convhull(Pre_X0),:);
p_Pre_X0 = polyshape(Pre_X0);
plot(p_Pre_X0);

% -------------------------------------------------------------------------

% print the C_inf values (looks right)
[C, d] = vert2con(Pre_X0);
Ans = [C, d];
Ans = Ans([1,6,4,5,2,3],:);
Ans = Ans .* [4;4;2.22;2.22;10;10]