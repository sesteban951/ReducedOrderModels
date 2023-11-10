%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fly wheel inverted pendulum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% system parameters
Ip = .5;
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
t_final = 20;
tspan = [0:.03: t_final];

% desired state
x_d = [0;     % theta pole
       0;     % theta flywheel
       0;     % velocity pole
       0;     % velocity flywheel
       0];    % aux input 

% define intital conditions
th_0 = pi;
th_dot_0 = 0;

x0 = [th_0;        % theta pole
      0;           % theta flywheel
      th_dot_0;    % velocity pole
      0;           % velocity flywheel
      0];          % aux input         

% control lyapunov function computaions
Q = diag([10 5 1])
R = 0.01;
F = [0 1 0; 0 0 1; 0 0 0];
G = [0; 0; 1];
% kp = 3;
% kd = 1;
% kdd = 5;
% K = [kp, kd, kdd];
% p = [-3,-1,-.5];
% K = place(F,G,p);
% P = lyap(Acl',Q);
[K, P, ~] = lqr(F,G,Q,R);
Acl = F - G*K;
eigs(Acl);
lam_max_P = max(eigs(P));
lam_min_P = min(eigs(P));
lam_max_Q = max(eigs(Q));
lam_min_Q = min(eigs(Q));

M = Q+2*P*G*K;
eigs(M)

% forward integrate
[t,x] = ode45(@(t,x)dynamics(t,x,x_d,params,P,K), tspan, x0); 

% compute extra stuff
delta_list = [];
V_list = []';
V_upper = [];
V_lower = [];
Vdot_list = [];
Vdot_upper = [];

for i = 1:length(t)
    
    % Torques
    x_ = x(i,:);
    th = x_(1);
    th_dot = x_(3);
    psi_dot = x_(4);
    v = x_(5);
    sig3 = m2*L^2 + m1*l^2 + Ip;
%     th_ddot= (g*l*m1*sin(alph+th)-bp*th_dot+L*g*m2*sin(th))/sig3...
%               + (bf/sig3)*psi_dot - v/sig3;
    th_ddot= (g*l*m1*sin(     th)-bp*th_dot+L*g*m2*sin(th))/sig3...
              + (bf/sig3)*psi_dot - v/sig3;

    [u,delta] = CLF_CBF_control(t(i),x_,x_d,params,P,K)
    delta_list = [delta_list, delta];

    % Lyapunov
    eta1 = th - x_d(1);
    eta2 = th_dot;
    eta3 = th_ddot;
    eta = [eta1; eta2; eta3];
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
          "Flywheel Vel., $\dot{\psi}$ [rad/s]";
          "Input, $v$ [Nm]"];

% plot all states
plt_rows = 3;
plt_cols = 3;
figure(1);
for i = 1:length(lables)
    subplot(plt_rows,plt_cols,i)
    plot(t, x(:,i),"Color",'blue','LineWidth',2)
    xlabel("Time, $t$ [sec]",Interpreter="latex")
    ylabel(lables(i),Interpreter="latex")
%     set(gcf, 'Position', [0 0 1200 1200])
    set(gca,'TickLabelInterpreter','latex')
    grid on;
    hold off;
end

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

tic
i = 1;
while i <= length(t)
    while toc<t(i)
        %
    end
    
    draw(t(i), x(i,:), params)
    i = i+1;
end

% core dynamics of flywheel inverted pendulum
function  xdot = dynamics(t,x,x_d,params,P,K)

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
    v = x(5)

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
    q_ddot = inv(D)*(-C*q_dot - G + B*v);

    % control affine form 
    f = [q_dot;
         inv(D) * (-C*q_dot - G)];
    g = [zeros(2,1);
         inv(D)*B];

    f_tild = [f + g*v;
              0];
    g_tild = [zeros(4,1);
              1];

    % CLF-CBF control
    [u, ~] = CLF_CBF_control(t,x,x_d,params,P,K);

    % return dynamics
    xdot = f_tild + g_tild*u;
    
end

%--------------------------------------------------------------------------

% CLF CBF QP controller
function [u, delta] = CLF_CBF_control(t,x,x_d,params,P,K)

    % states
    th = x(1);
    psi = x(2);         % position of the flywheel doesn't affect dynamics
    th_dot = x(3);
    psi_dot = x(4);
    v = x(5);

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

%     th_ddot= (g*l*m1*sin(alph+th)-bp*th_dot+L*g*m2*sin(th))/sig3...
%               + (bf/sig3)*psi_dot - v/sig3;
    th_ddot= (g*l*m1*sin(     th)-bp*th_dot+L*g*m2*sin(th))/sig3...
              + (bf/sig3)*psi_dot - v/sig3;

%     psi_ddot = (-g*l*m1*sin(alph+th)+bp*th_dot-L*g*m2*sin(th))/sig3...
%                -(bf*sig4/sig1)*psi_dot + (sig4/sig1) *v;
    psi_ddot = (-g*l*m1*sin(     th)+bp*th_dot-L*g*m2*sin(th))/sig3...
               -(bf*sig4/sig1)*psi_dot + (sig4/sig1) *v;

    % raw output dynamic
    f_hat_1 = th_dot;
    f_hat_2 = th_ddot;
%     f_hat_3 = (g*l*m1*cos(alph+th)*th_dot-bp*th_ddot+L*g*m2*cos(th)*th_dot)/sig3...
%               + (bf/sig3)*psi_ddot;
    f_hat_3 = (g*l*m1*cos(     th)*th_dot-bp*th_ddot+L*g*m2*cos(th)*th_dot)/sig3...
              + (bf/sig3)*psi_ddot;
    f_hat = [f_hat_1; f_hat_2; f_hat_3];

    g_hat_1 = 0;
    g_hat_2 = 0;
    g_hat_3 = -1/sig3;
    g_hat = [g_hat_1; g_hat_2; g_hat_3];

    % define output dynmics coordinates
    eta1 = th - x_d(1);
    eta2 = th_dot;
    eta3 = th_ddot;
    eta = [eta1; eta2; eta3]

    % Lyapunov
    lam = 1;
    V = eta'*P*eta;
%     LfV = f_hat'*P*eta + eta'*P*f_hat
%     LgV = 2*eta'*P*g_hat
    F = [0 1 0; 0 0 1; 0 0 0];
    G = [0; 0; 1];
    LFV = eta'*(F'*P + P*F)*eta;
    LGV = 2*eta'*P*G;

    % CLF QP
%     A = LGV;
%     b = -lam*V - LFV;
%     H = 1;
%     F = 1;
%     t
%     w = quadprog(H,F,A,b)
% %     w = -K*eta;
%     
%     u = pinv(g_hat(3))*(-f_hat(3) + w);
%     delta =0;

%     %%%%%%%%%%%%%%%%%%%%%%%%%%% CBF Part %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define control barrier functions
    beta = 10;

    w_max = 20;    % Eric max torque is 2.1 Nm pretty low
    w_min = -20;
    h1 = w_max - v;
    h2 = v - w_min;

    Lfh1 = 0;
    Lgh1 =-1;

    Lfh2 = 0;
    Lgh2 = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%% QP  Part %%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_scale = 10;
    delta_scale = 10000;
    H = [u_scale, 0;
         0, delta_scale];
    F = [u_scale;
         delta_scale];

    A = [LGV,  -1;
        -Lgh1, 0;
        -Lgh2, 0];
    b = [-lam*V - LFV;
         beta*h1 + Lfh1;
         beta*h2 + Lfh2];

    sol = quadprog(H,F,A,b);

    w = sol(1);
    delta = sol(2);

    t
    u = pinv(g_hat(3))*(-f_hat(3) + w);

end
