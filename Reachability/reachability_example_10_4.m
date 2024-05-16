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
U = [5; -5];

% ---------------------- Forwards Reachability -----------------------------
figure(); grid on; hold on;

% plot the orignal intial condition set
X0 = con2vert(H_x0, h_x0);
X0 = X0(convhull(X0),:);
p_X0 = polyshape(X0);
plot(p_X0);

% plot the one step reachable set
AX = X0 * A';
[AX_row, ~] = size(AX);
BU = U * B';
[BU_row, ~] = size(BU);

AX_sum_BU = [];
for i = 1:AX_row
    for j = 1:BU_row
        AX_sum_BU = [AX_sum_BU; AX(i,:) + BU(j,:)]
    end
end

AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);

[C, d] = vert2con(AX_sum_BU)
[row, ~] = size(C);
C = [C; H_x];
d = [d; h_x];
AX_sum_BU = con2vert(C, d);
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% plot another one!
AX =  AX_sum_BU * A';
[AX_row, ~] = size(AX);
BU = U * B';
[BU_row, ~] = size(BU);

AX_sum_BU = [];
for i = 1:AX_row
    for j = 1:BU_row
        AX_sum_BU = [AX_sum_BU; AX(i,:) + BU(j,:)]
    end
end

AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
[C, d] = vert2con(AX_sum_BU)
[row, ~] = size(C);
C = [C; H_x];
d = [d; h_x];
AX_sum_BU = con2vert(C, d);
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% plot another one!
AX =  AX_sum_BU * A';
[AX_row, ~] = size(AX);
BU = U * B';
[BU_row, ~] = size(BU);

AX_sum_BU = [];
for i = 1:AX_row
    for j = 1:BU_row
        AX_sum_BU = [AX_sum_BU; AX(i,:) + BU(j,:)]
    end
end

AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
[C, d] = vert2con(AX_sum_BU)
[row, ~] = size(C);
C = [C; H_x];
d = [d; h_x];
AX_sum_BU = con2vert(C, d);
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% plot another one!
AX =  AX_sum_BU * A';
[AX_row, ~] = size(AX);
BU = U * B';
[BU_row, ~] = size(BU);

AX_sum_BU = [];
for i = 1:AX_row
    for j = 1:BU_row
        AX_sum_BU = [AX_sum_BU; AX(i,:) + BU(j,:)]
    end
end

AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
[C, d] = vert2con(AX_sum_BU)
[row, ~] = size(C);
C = [C; H_x];
d = [d; h_x];
AX_sum_BU = con2vert(C, d);
AX_sum_BU = AX_sum_BU(convhull(AX_sum_BU),:);
p_AX_sum_BU = polyshape(AX_sum_BU);
plot(p_AX_sum_BU);

% --------------------------------------------------------------------
legend(["$\mathcal X_0$", "$\mathcal R_1 (\mathcal X_0)$", "$\mathcal R_2 (\mathcal X_0)$", "$\mathcal R_3 (\mathcal X_0)$", "$\mathcal R_4 (\mathcal X_0)$"], 'interpreter', 'latex');