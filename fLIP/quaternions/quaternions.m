clear;clc;
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
plot_quat(quat_1,'ro',[a;b;c]);

%% error
global error_cross
error_cross = @(q) [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];
otimes = @(q) [q(1) -q(2) -q(3) -q(4); 
                q(2) q(1) -q(4) q(3);
                q(3) q(4) q(1) -q(2);
                q(4) -q(3) q(2) q(1)];
Exp = @(w) [cos(norm(w)/2); w/norm(w)*sin(norm(w)/2)];    
Log = @(q) 2*q(2:4)*atan2(norm(q(2:4)),q(1))/norm(q(2:4));


d_n = q_0(1)*q_1(1)+q_0(2:4)*q_1(2:4)';
d_q = q_0(1)*q_1(2:4)' - q_1(1)*q_0(2:4)'-error_cross(q_0(2:4))*q_1(2:4)';

%% dynamics
clf
omega = [0 1 0]';
view(50,30)
t_max = 2;
[t,q] = ode45(@(t,x) quatODE(t,x),[0 t_max],q_0');
q_old = q;
t_old = t;
q_dot_old = [diff(q_old,1)./diff(t_old); 0 0 0 0];

t = linspace(0,t_max,20);
dt = mean(diff(t));
q = interp1(t_old,q_old,t);
q_dot = interp1(t_old,q_dot_old,t);
omega = zeros(length(t),3);
q_dot_pred = zeros(length(t),4);
for i = 1:length(t)
    [~,o] = quatODE(t(i),q(i,:)');
    zeta_k(i,:) = Log(q(i,:)); % The MPC decision variables
    q_k(i,:) = Exp(zeta_k(i,:)');
    omega(i,:) = o;
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
plot(t,q_dot_pred,'--')
title('\dot{Quat} vs \dot{Quat} pred')
subplot(2,3,4)
plot(t,zeta_k)
hold on;
plot(t+dt,zeta_kp1,'o--')
title('Log(Quat) vs Log(Quat) pred')
subplot(2,3,[5 6])
while true
    tmp = toc;
    ind = find(t>tmp,1);
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
function [dq, omega] = quatODE(t,q)
global error_cross
K = diag([5 5 5]);
% q_d = [sqrt(2)/2 0 sqrt(2)/2 0];
q_d = [1 0 0 0];

delta_q = q(1)*q_d(2:4)' - q_d(1)*q(2:4)-error_cross(q(2:4))*q_d(2:4)';
e_o = delta_q;
omega = K*e_o;
dq = 1/2*[0 -omega'; omega -error_cross(omega)]*q;
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