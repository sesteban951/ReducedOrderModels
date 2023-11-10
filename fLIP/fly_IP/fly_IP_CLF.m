%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fly wheel inverted pendulum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; clear figure;

% system parameters
Ip = 1;   % hopper torso inertia
If = 1;  % flywheel moment of inertia
bp = 0;
bf = 0.01;   % 
m1 = 1.0; % torso mass
m2 = 0.5; % flywheel mass
L = 1; % hopper height
l = L;     % CG off set along pole
alph = 0.5;      % CG off set along pole
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
tspan = [0:.03:t_final];

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
kp = 1;
kd = 1;
K = [kp kd];         % fbl closed loop gains
Acl = (F - G*K);   % closed loop fbl system
Q = eye(eta_dim);  % CTLE
P = lyap(Acl',Q);  % CTLE, you have to give lyap A^T

% Vdot <= 0 stable <==> (Q+2PGK) is P.S.D
M = Q + 2*P*G*K
eigs_M = eigs(M)
M_PSD = all(eigs_M >=  0)

lam_min_P = min(eig(P));  % real symmetric matrices have real eigenvalues
lam_max_P = max(eig(P));
lam_min_Q = min(eig(Q));
lams = [lam_min_P, lam_max_P, lam_min_Q];

% forward integrate
[t,x] = ode45(@(t,x)dynamics(t,x,x_d,params,F,G,P,K,lams), tspan, x0); 

% compute extra stuff
u_list = [];
V_list = []';
V_upper = [];
V_lower = [];
Vdot_list = [];
Vdot_upper = [];
for i = 1:length(t)
    
    % Torques
    x_ = x(i,:);
    u = CLF_control(t(i),x_,x_d,params,F,G,P,K,lams);
    u_list = [u_list, u];

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
lables = ["Angle, $\theta$ [rad]";
          "Angle, $\psi$ [rad]";
          "Vel., $\dot{\theta}$ [rad/s]";
          "Vel., $\dot{\psi}$ [rad/s]"];

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
legend(["$V$","$V_{lower}$","$V_{upper}$"], Interpreter='latex')
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
legend(["$V$","$\dot{V}_{upper}$"], Interpreter='latex')
% set(gcf, 'Position', [0 0 1200 1200])
set(gca,'TickLabelInterpreter','latex')
grid on;

%% 
% draw fLIP
figure;
[row, col] = size(x);
tic
i=1;
while i <= length(t)
    while toc < t(i)
        %
    end
    draw(t(i),x(i,:),params)
    toc;
    t(i)
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUXILLARY 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%     G = [-m1*g*l*sin(th     ) - m2*g*L*sin(th);
%          0];
    B = [0;
         1];

    % configuration
    q_dot = [th_dot; 
             psi_dot];

    % CLF control
    u = CLF_control(t,x,x_d,params,F,G,P,K,lams);

    % control affine form 
    f = [q_dot;
         inv(D) * (-C*q_dot - G)];
    g = [zeros(2,1);
         inv(D)*B];

    % return dynamics
    xdot = f + g*u;
    
end

% CLF controller
function u = CLF_control(~,x,x_d,params,F,G,P,K,lams)

    % states
    th = x(1);
    psi = x(2);          % position of the flywheel doesnt affect dynamics
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

    % define output dynmics coordinates
    eta1 = th - x_d(1);
    eta2 = th_dot;
    eta = [eta1; eta2];

    % controller
    v = -K*eta;
    fbl_part = -(g*l*m1*sin(alph+th)-bp*th_dot+L*g*m2*sin(th))/sig3...
               - (bf/sig3)*psi_dot;  % conundrum: we dont acutall know alpha
%     fbl_part = -(g*l*m1*sin(      th)-bp*th_dot+L*g*m2*sin(th))/sig3...
%                - (bf/sig3)*psi_dot; 
    u = -sig3*( fbl_part + v);

end

%--------------------------------------------------------------------------

% core dynamics of flywheel inverted pendulum, W/ INTEGRAL STATE
function  xdot = dynamics_e(t,x,x_d,params)

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
    err1 = th - x_d(5);
    err2 = th_dot - x_d(6);

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

    % PID control
    u = PID_control(t,x, x_d, err1, err2);

    % control affine form 
    f = [q_dot;
         inv(D) * (-C*q_dot - G)];
    g = [zeros(2,1);
         inv(D)*B];

    % return dynamics
    xdot = zeros(6,1);
    xdot(1:4) = f + g*u;
    xdot(5:6) = [err1; err2];
    
end

% PID cotnrol
function u = PID_control(t,x,x_d, err1, err2)
    
    % states
    th = x(1);
    psi = x(2);              % position of the flywheel doesnt affect dynamics
    th_dot = x(3);
    psi_dot = x(4);

    % states desired
    th_d = x_d(1);
    psi_d = x_d(2);          % position of the flywheel doesnt affect dynamics
    th_dot_d = x_d(3);
    psi_dot_d = x_d(4);

    % PID gains
    P = 60;
    I = 10;
    D = 6;

    % if close to eq point use PID
    err_norm = norm([err1; err2])
    if err_norm < .1
        u = P*err1 ...
            + I*err1 ...
            + D*err2 ;
        fprintf("using PID\n")
    
    % if not close to the equilibrium point use PD
    else
        u = P*(err1) + D*(err2) ;
        fprintf("using PD\n")
        
    end
    
    % add staturation
    u_max = 50;
    u_min = -50;

    if u > u_max
        u = u_max;
    elseif u < u_min
        u = u_min;
    end

end





