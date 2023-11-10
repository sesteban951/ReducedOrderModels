clear;clc;
init();
a = [1,0,0];
b = [0,1,0];
c = [0,0,1];
% quat = quaternion([30,45,0],'eulerd','ZYX','point');
% q_0 = [1 0 0 0];
q_0 = randrot().compact;
q_1 = [sqrt(2)/2 0 sqrt(2)/2 0];

quat_0 = quaternion(q_0);
quat_1 = quaternion(q_1);

[p,l] = plot_quat(quat_0,'bo',[a;b;c]);

%% error
global error_cross Exp
error_cross = @(q) [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];
otimes = @(q) [q(1) -q(2) -q(3) -q(4); 
                q(2) q(1) -q(4) q(3);
                q(3) q(4) q(1) -q(2);
                q(4) -q(3) q(2) q(1)];
Exp = @(w) [cos(norm(w)/2); w/norm(w)*sin(norm(w)/2)];    
Log = @(q) 2*q(2:4)*atan2(norm(q(2:4)),q(1))/norm(q(2:4));

%% MPC
N = 50;
x0 = [q_0 0 0 0];
X0 = [Log(q_0) 0 0 0]';
Xf = [0 0 0 0 0 0]';
t_max = 5;
tspan = [0 t_max];
dt = 0.1;
u_min = -200;
u_max = 200;
state_stage_cost = 10;
input_stage_cost = 1;

%%% MPC size
state_dim = 6;
state_size_without_terminal = state_dim*(N-1);
state_size_with_terminal = state_dim*N;
state_end_index=state_dim*N;
input_dim = 3;
input_size = input_dim*(N-1);
input_end_index=state_end_index+input_dim*(N-1);
total_dim = state_dim*N + input_dim*(N-1);

% Pre-allocate program
X_MPC = [];
T_MPC = [];
U_MPC = [];
X_K_MPC_CLF = [];
U_FF_MPC_CLF = [];
Ad = [eye(3) dt*eye(3); zeros(3) eye(3)];
Bd = [zeros(3); eye(3)*dt];

% Optimizer variables
states = sdpvar(state_dim,N);
inputs = sdpvar(input_dim,N-1);
x_lin_ = sdpvar(N-1,state_dim);
u_lin_ = sdpvar(N-1,input_dim);
x0_ = sdpvar(state_dim,1);

Constraints = [];
% Decision variables
Constraints = [Constraints states(:,1) == x0_];
Constraints = [Constraints states(:,end) == Xf];

for i = 1:N-1
    % Dynamics
    Constraints = [Constraints states(:,i+1) == Ad*states(:,i) + Bd*inputs(:,i)];
end

Q = state_stage_cost*eye(state_dim*N);
R = input_stage_cost*eye(input_dim*(N-1));

Objective = 1/2*(states(:)'*Q*states(:) + inputs(:)'*R*inputs(:));

P = optimizer(Constraints,Objective,sdpsettings('solver','mosek','verbose',0),...
    {x0_},{states,inputs});

x = x0;

for iter = 1:tspan(end)/dt
    disp(iter)
    x0 = x(end,:)';
    frak_x0 = [Log(x0(1:4)); x0(5:7)];
    [sol, diagnostics,d1,d2,d3,d4] = P({frak_x0});
    if diagnostics ~= 0
            error('Issue with Mosek in proposed');
    end
    
    t_MPC = 0:dt:dt*(N-1);
    x_MPC = sol{1}';
    u_MPC = sol{2}';
    
    [t,x] = ode45(@(t,x) quatODE(t,x,x_MPC(2,:), u_MPC(2,:)),[0 dt],x0');
    
    if isempty(T_MPC)
        T_MPC = [T_MPC; t];
        X_MPC = [X_MPC; x];
        U_MPC = [U_MPC; repmat(u_MPC(2,:),length(t)-1,1)];
    else
        T_MPC = [T_MPC; t(2:end)+T_MPC(end)];
        X_MPC = [X_MPC; x(2:end,:)];
        U_MPC = [U_MPC; repmat(u_MPC(2,:),length(t)-1,1)];
    end
    
end
% plot(T_MPC,X_MPC)
% pause();


%% dynamics
clf
view(50,30)

x0 = [q_0 0 0 0];

t = T_MPC;
x = X_MPC;

q = x(:,1:4);
w = x(:,5:7);

x_old = x;
q_old = q;
w_old = w;
t_old = t;
q_dot_old = [diff(q_old,1)./diff(t_old); 0 0 0 0];

t = linspace(0,t_max,t_max/dt);
x = interp1(t_old,x_old,t);
q = interp1(t_old,q_old,t);
q_dot = interp1(t_old,q_dot_old,t);
omega = interp1(t_old,w_old,t);
q_dot_pred = zeros(length(t),4);
for i = 1:length(t)
%     [~,acc] = quatODE(t(i),x(i,:)');
%     tau(i,:) = acc';
    zeta_k(i,:) = Log(q(i,:)); % The MPC decision variables
    q_k(i,:) = Exp(zeta_k(i,:)');
    q_dot_pred(i,:) = (1/2*otimes(q_k(i,:))*[0;omega(i,:)']);
    q_pred(i,:) = otimes(q_k(i,:))*Exp(omega(i,:)'*dt);
    zeta_kp1(i,:) = zeta_k(i,:) + omega(i,:)*dt;
end



tic

subplot(2,3,1)
plot(t,q)
hold on;
plot(t+dt,q_pred,'o--')
title('Quat vs Quat pred')
subplot(2,3,2)
plot(t,omega)
title('Omega (PD feedback)')
subplot(2,3,3)
plot(t,q_dot)
hold on;
plot(t,q_dot_pred,'o--')
title('\dot{Quat} vs \dot{Quat} pred')
subplot(2,3,4)
plot(t,zeta_k)
hold on;
plot(t+dt,zeta_kp1,'o--')
title('Log(Quat) vs Log(Quat) pred')
subplot(2,3,[5 6])
while true
    tmp = toc;
    ind = find(t(1:end-1)>tmp,1);
    if isempty(ind)
        break;
    else
        view(50,30)
%         delete(p);
        delete(l);
    end
    [p,l] = plot_quat(quaternion(q(ind,:)),'bo',[a;b;c]);
    drawnow
end

%%
function [dx, tau] = quatODE(t,x, x_MPC, u_MPC)
global error_cross Exp

q = x(1:4);
w = x(5:7);

K_p = diag([50 50 50]);
K_d = diag([10 10 10]);
% q_d = [sqrt(2)/2 0 sqrt(2)/2 0];
% q_d = [1 0 0 0];
q_d = Exp(x_MPC(1:3)')';
w_d = x_MPC(4:6)';

delta_q = q(1)*q_d(2:4)' - q_d(1)*q(2:4)-error_cross(q(2:4))*q_d(2:4)';
e_o = delta_q;
tau = K_p*e_o - K_d*(w-w_d);

dq = 1/2*[0 -w'; w -error_cross(w)]*q;
dw = tau+u_MPC';

dx = [dq; dw];
end

function [p,l] = plot_quat(quat,color,frame)

rP = rotateframe(quat,frame);
pos = [0 0 0];

hold on;
p(1) = plot3(rP(1,1),rP(1,2),rP(1,3),color);
p(2) = plot3(rP(2,1),rP(2,2),rP(2,3),color);
p(3) = plot3(rP(3,1),rP(3,2),rP(3,3),color);
l(1) = plot3([0;rP(1,1)],[0;rP(1,2)],[0;rP(1,3)],'k');
l(2) = plot3([0;rP(2,1)],[0;rP(2,2)],[0;rP(2,3)],'k');
l(3) = plot3([0;rP(3,1)],[0;rP(3,2)],[0;rP(3,3)],'k');
% hold off

grid on
axis([-1 1 -1 1 -1 1])
xlabel('x')
ylabel('y')
zlabel('z')
end