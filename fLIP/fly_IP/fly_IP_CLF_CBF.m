%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fly wheel inverted pendulum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% system parameters
Ip = 0.5;
If = 2;
bp = 0.1;
bf = 0.1;
m1 = .5;
m2 = .25;
L = 1.5;
l = .9*L;     % CG off set along pole
alph = 0.1;      % CG off set along pole
g = 9.81;

% store system parameters 
params = [Ip;
          If;
          bp;
          bf;
          m1;
          m2;
          L;
          l;
          alph;
          g];

% define time span
t_final = 25;
tspan = [0 t_final];

% desired state
x_d = [0;     % theta pole
       0;     % theta flywheel
       0;     % velocity pole
       0];     % velocity flywheel
%        0;     % theta error
%        0];    % pole velocity error

% define intital conditions

th_0 = pi;
th_dot_0 = 0;

x0 = [th_0;        % theta pole
      0;           % theta flywheel
      th_dot_0;    % velocity pole
      0];           % velocity flywheel
%       th_0 - x_d();        % theta error
%       th_dot_0];   % pole velocity error

% control lyapunov function computaions
eta_dim = 2;
F = [0 1; 0 0];    % feedback linearization
G = [0 ; 1];       % feedback linearization
kp = 5;
kd = 5;
K = [kp kd];         % fbl closed loop gains
Acl = (F - G*K);   % closed loop fbl system
Q = eye(eta_dim);  % CTLE
P = lyap(Acl',Q);  % CTLE, you have to give lyap A^T

% Vdot <= 0 stable <==> (Q+2PGK) is P.S.D.
M = Q + 2*P*G*K
eigs_M = eigs(M)
M_PSD = all(eigs_M >=  0) % P.S.D. by checking is all eigs >= 0

lam_min_P = min(eig(P));  % real symmetric matrices have real eigenvalues
lam_max_P = max(eig(P));
lam_min_Q = min(eig(Q));
lams = [lam_min_P, lam_max_P, lam_min_Q];

% forward integrate
[t,x] = ode45(@(t,x)dynamics(t,x,x_d,params,F,G,P,K,lams), tspan, x0); 

% compute extra stuff
u_list = [];
delta_list = [];
V_list = []';
V_upper = [];
V_lower = [];
Vdot_list = [];
Vdot_upper = [];
for i = 1:length(t)
    
    % Torques
    x_ = x(i,:);
    [u,delta] = CLF_CBF_control(t(i),x_,x_d,params,F,G,P,K,lams);
    u_list = [u_list, u];
    delta_list = [delta_list, delta];

    % Lyapunov
    eta1 = x_(1) - x_d(1);
    eta2 = x_(3);
    eta = [eta1; eta2];
    V = eta'* P* eta;
    V_list = [V_list, V];
    V_up = lam_max_P*norm(eta)^2;
    V_upper = [V_upper, V_up];
    V_low = lam_min_P*norm(eta)^2;
    V_lower = [V_lower, V_low];

    Vdot = -eta' * M * eta;
    Vdot_list = [Vdot_list, Vdot];
    Vdot_up = - (lam_min_Q/lam_max_P) * V;
    Vdot_upper = [Vdot_upper, Vdot_up];

end

% plot results
lables = ["Pole Angle, $\theta$ [rad]";
          "Flywheel Angle, $\psi$ [rad]";
          "Pole Vel., $\dot{\theta}$ [rad/s]";
          "Flywheel Vel., $\dot{\psi}$ [rad/s]"];

% plot all states
plt_rows = 3;
plt_cols = 3;
figure(1);
for i = 1:4
    subplot(plt_rows,plt_cols,i)
    plot(t, x(:,i),"Color",'blue','LineWidth',2)
    xlabel("Time, $t$ [sec]",Interpreter="latex")
    ylabel(lables(i),Interpreter="latex")
%     set(gcf, 'Position', [0 0 1200 1200])
    set(gca,'TickLabelInterpreter','latex')
    grid on;
    hold off;
end

% plot inputs
i = i+1;
subplot(plt_rows,plt_cols, i)
plot(t, u_list,"Color",'blue','LineWidth',2)
xlabel("Time, $t$ [sec]",Interpreter="latex")
ylabel("Input, $u$",Interpreter="latex")
set(gca,'TickLabelInterpreter','latex')
grid on;

% plot relaxation
i = i+1;
subplot(plt_rows,plt_cols, i)
plot(t, delta_list,"Color",'blue','LineWidth',2)
xlabel("Time, $t$ [sec]",Interpreter="latex")
ylabel("Relaxation, $\delta$",Interpreter="latex")
set(gca,'TickLabelInterpreter','latex')
grid on;

% plot V
i = i+1;
subplot(plt_rows,plt_cols,i)
hold on;
plot(t,V_list,"Color","black",'LineWidth',2)
plot(t,V_lower,"Color","blue",'LineWidth',2)
plot(t,V_upper,"Color","red",'LineWidth',2)
hold off;
xlabel("Time, $t$ [sec]",Interpreter="latex")
ylabel("$V$",Interpreter="latex")
% set(gcf, 'Position', [0 0 1200 1200])
set(gca,'TickLabelInterpreter','latex')
grid on;

% plot Vdot
i = i+1;
subplot(plt_rows,plt_cols,i)
hold on;
plot(t,Vdot_list,"Color","black",'LineWidth',2)
plot(t,Vdot_upper,"Color","red",'LineWidth',2)
hold off;
xlabel("Time, $t$ [sec]",Interpreter="latex")
ylabel("$\dot{V}$",Interpreter="latex")
% set(gcf, 'Position', [0 0 1200 1200])
set(gca,'TickLabelInterpreter','latex')
grid on;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw fLIP
figure(2); 
[row, col] = size(x);

size_des = 600;     % has to be less than total data size
n = size(t)/size_des;
n = round(n);
T_ = t(1:n:end);
X_ = x(1:n:end,:);

tic;
i = 1;
while i < length(T_)
    while toc < T_(i)
        %
    end
    draw(T_(i), X_(i,:), params)
    i = i+1;
end

% core dynamics of flywheel inverted pendulum
function  xdot = dynamics(t,x,x_d,params,F,G,P,K,lams)

    % unpack system parameters
    Ip = params(1);
    If = params(2);
    bp = params(3);
    bf = params(4);
    m1 = params(5);
    m2 = params(6);
    L = params(7);
    l = params(8);       % CG off set along pole
    alph = params(9);    % CG off set along pole
    g = params(10);

    % states
    th = x(1);
    psi = x(2);           % position of the flywheel doesnt affect dynamics
    th_dot = x(3);
    psi_dot = x(4);

    % Mech Sys. Dynamics 
    D = [Ip + m1*l^2 + If + m2*L^2, If;
         If,                        If];
    C = [bp, 0;
         0 , bf];
    G = [-m1*g*l*sin(th+alph) - m2*g*L*sin(th);
         0];
    B = [0;
         1];

    % configuration
    q_dot = [th_dot; 
             psi_dot];

    % CLF-CBF control
    [u, ~] = CLF_CBF_control(t,x,x_d,params,F,G,P,K,lams);

    % control affine form 
    f = [q_dot;
         inv(D) * (-C*q_dot - G)];
    g = [zeros(2,1);
         inv(D)*B];

    % return dynamics
    xdot = f + g*u;
    
end

%--------------------------------------------------------------------------

% CLF CBF QP controller
function [u, delta] = CLF_CBF_control(t,x,x_d,params,F,G,P,K,lams)

    % states
    th = x(1);
    psi = x(2);          % position of the flywheel doesn't affect dynamics
    th_dot = x(3);
    psi_dot = x(4);

    % system parameters
    Ip = params(1);
    If = params(2);
    bp = params(3);
    bf = params(4);
    m1 = params(5);
    m2 = params(6);
    L = params(7);
    l = params(8);
    alph = params(9);
    g = params(10);

    sig1 = If*m2*L^2 + If*m1*l^2 + If*Ip;
    sig3 = m2*L^2 + m1*l^2 + Ip;
    sig4 = m2*L^2 + m1*l^2 + Ip + If;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% CLF Part %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % raw output dynamic
    f_hat_1 = th_dot;
    f_hat_2 = (g*l*m1*sin(alph+th)-bp*th_dot+L*g*m2*sin(th))/sig3...
              + (bf/sig3)*psi_dot;
    f_hat_2 = (g*l*m1*sin(     th)-bp*th_dot+L*g*m2*sin(th))/sig3...
              + (bf/sig3)*psi_dot; % dont actuall yknow alpha
    f_hat = [f_hat_1; f_hat_2];

    g_hat_1 = 0;
    g_hat_2 = -1/sig3;
    g_hat = [g_hat_1; g_hat_2];

    % define output dynmics coordinates
    eta1 = th - x_d(1);
    eta2 = th_dot;
    eta = [eta1; eta2];

    % Lyapunov
    lam = 1;
    V = eta'*P*eta;
    LfV = f_hat'*P*eta + eta'*P*f_hat;
    LgV = 2*eta'*P*g_hat;

    % CLF QP
%     A = LgV;
%     b = -lam*V - LfV;
%     H = 1;
%     F = 1;
% 
%     t
%     u = quadprog(H,F,A,b)

%     %%%%%%%%%%%%%%%%%%%%%%%%%%% CBF Part %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define control barrier functions
    beta = 1;

    w_max = 17;    % right before it drops off in eric's diagram
    w_min = -17;
    h1 = w_max - psi_dot;
    h2 = psi_dot - w_min;

    Lfh2 = (-g*l*m1*sin(alph+th) + bp*th_dot - L*g*m2*sin(th))/sig3 ...
           - (bf*sig4/sig1)*psi_dot;
    Lfh2 = (-g*l*m1*sin(     th) + bp*th_dot - L*g*m2*sin(th))/sig3 ...
           - (bf*sig4/sig1)*psi_dot; % dont actually know alpha
    Lgh2 = sig4/sig1;

    Lfh1 = -Lfh2;
    Lgh1 = -Lgh2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% QP  Part %%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_scale = 10;
    delta_scale = 10000;
    H = [u_scale, 0;
         0, delta_scale];
    F = [u_scale;
         delta_scale];

    A = [LgV,  -1;
        -Lgh1, 0;
        -Lgh2, 0];
    b = [-lam*V - LfV;
         beta*h1 + Lfh1;
         beta*h2 + Lfh2];

    sol = quadprog(H,F,A,b);

    u = sol(1);
    delta = sol(2);

end
