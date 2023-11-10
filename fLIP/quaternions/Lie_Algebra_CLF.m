% testing CLF on a Lie algebra
clear all; clf; clc;

% inital and final conditions
% q_0 = randrot().compact; % generate random unit quaternion
% q_0 = q_0'
q_0 = [-.5 -.1 -.4 -.2]';
q_0 = q_0 / norm(q_0);
% q_0 = [-0.7384; 0.6243; 0.1939; -0.1656];
% w_0 = 2 * ([rand(1), rand, 0]' - 0.5) % genrate random angular speed
w_0 = [0;0;0];

% inital/final conditions in euclidean space
x_0 = [q_0; w_0];
% tspan = [0, 20];
tspan = [0:.03:20.5];

% go to the identity element with zero velocity
q_d = [0 0 1 0]';  % desired postion 
w_d = [0 0 0]';    % desired velocity
x_d = [q_d; w_d];

% non inertia output dyanmics, V(eta) = 1/2 Eta' * P * eta
A = [zeros(3), eye(3); zeros(3), zeros(3)];
B = [zeros(3); eye(3)];
% C = ctrb(A,B);
% C_rank = rank(C);
% Qx = diag([10, 10, 10, 1, 1, 1]);
% Qu = diag([1, 1, 1]);
% K = lqr(A,B,Qx,Qu);
kp = 1.5;    % can directly specify PD gains
kd = .5;    % can directly specify PD gains
K = [diag(kp*[1 1 1]) diag(kd*[1 1 1])]; % can directly specify PD gains
% poles = [-6 -5 -4 -3 -2 -1];
% [K, ~] = place(A, B, poles)
Acl = A - B*K;
Acl_eig = eig(Acl)
Q = 1*diag([1 1 1 1 1 1]);
P = lyap(Acl',Q);  % you have to give lyap A^T
% P = icare(A,B, Q);
% CTLE_lyap = Acl'*P_lyap + P_lyap*Acl;
% CTLE_care = A'*P_care + P_care*A - P_care*B*B'*P_care;

% compute lyaupnov lambdas
lam_min_P = min(eig(P));
lam_max_P = max(eig(P));
lam_min_Q = min(eig(Q));
lams = [lam_min_P, lam_max_P, lam_min_Q];

% forward integrate dynamics
% [t, x] = ode45(@(t,x)quatODE_inertia(t,x,x_d), tspan, x_0);
% [t, x] = ode45(@(t,x)quatODE_PD(t,x,x_d), tspan, x_0);
[t, x] = ode45(@(t,x)quatODE_clf_QP(t,x,x_d, A,B,Q,P, lams), tspan, x_0);

% for plotting lyapunov function
quats = x(:,1:4);
w = x(:,5:7);

V_vals = [];
V_lowers = [];
V_uppers = [];

Vdot_vals = [];
Vdot_uppers = [];

U_vals = [];

for i = 1:length(quats)

    % compute Lyapunov vars
    q_t = quats(i,:)';
    e_pos = mult(q_conj(q_d), q_t);
    e_pos = mult(q_t, q_conj(q_d));
    xi_t = Log(e_pos);
    w_t = w(i,:)';
    eta_t = [xi_t;w_t];

    % compute V(xi,w)
    V = eta_t' * P * eta_t;
%     V = 0.5 * (xi' * xi) +  0.5 * (w' * w);
    V_lower = lam_min_P * norm(eta_t)^2;
    V_upper = lam_max_P * norm(eta_t)^2;
    V_lowers = [V_lowers V_lower];
    V_uppers = [V_uppers V_upper];
    V_vals = [V_vals V];

    % compute input
    u_t = -0.5*B'*P*eta_t;    % LQR solution
    U_vals = [U_vals; u_t'];

    Vdot = eta_t' * (-Q) * eta_t + 2 * eta_t' * P * B * u_t;
    Vdot_vals = [Vdot_vals, Vdot];
    Vdot_upper = - (lam_min_Q / lam_max_P) * V;
    Vdot_uppers = [Vdot_uppers, Vdot_upper];

end

quats(end)
w(end)

% plot quaternions states
figure(1)
subplot(2,3,1)
plot(t,x(:,1:4))
xlabel("Time, $t$", Interpreter="latex")
ylabel("Quaternion Component, $q_i$", Interpreter="latex")
legend(["$q_0$", "$q_1$", "$q_2$", "$q_3$"], Interpreter="Latex")
grid on;

% plot angular velocity states
subplot(2,3,2)
plot(t,x(:,5:7))
xlabel("Time, $t$", Interpreter="latex")
ylabel("Velocity Component, $\omega_i$", Interpreter="latex")
legend(["$w_1$", "$w_2$", "$w_3$"], Interpreter="Latex")
grid on;

% plot u
subplot(2,3,3)
plot(t,U_vals)
xlabel("Time, $t$", Interpreter="latex")
ylabel("$k(x)$", Interpreter="latex")
legend(["$u_1$", "$u_2$", "$u_3$"], Interpreter="Latex")
grid on;

% plot V(eta)
subplot(2,3,4)
hold on;
plot(t,V_lowers,'b')
plot(t,V_vals,'k')
plot(t,V_uppers,'r')
xlabel("Time, $t$", Interpreter="latex")
ylabel("$V$", Interpreter="latex")
legend(["$\lambda_{min}(P) \|\eta\|^2$", "$V(\eta)$","$\lambda_{max}(P) \|\eta\|^2$"], Interpreter="Latex")
grid on;

% plot V_dot(eta)
subplot(2,3,5)
hold on;
plot(t,Vdot_vals,'k')
plot(t,Vdot_uppers,'r')
xlabel("Time, $t$", Interpreter="latex")
ylabel("$\dot{V}$", Interpreter="latex")
legend(["$\dot{V}(\eta)$", "$-\frac{\lambda_{min}(Q)}{\lambda_{max}(P)} V(\eta)$"], Interpreter="Latex")
grid on;

% play animation of attitude control
pause(0.5);

tic
i=1;
while i<length(t)
    while toc < t(i)
        %
    end
    figure(2);
    quat = quaternion(quats(i,1), quats(i,2), quats(i,3), quats(i,4));
    str = sprintf("Time = %.3f",t(i));
    poseplot(quat);
    title(str);

    i =i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILLARY Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotating body dynamics
function xdot = quatODE_clf_QP(~, x, x_d, A,B,Q,P, lams)
    % get lambdas for CLF
    lam_max_P = lams(2);
    lam_min_Q = lams(3);

    % desired eq points
    q_d = x_d(1:4);
    w_d = x_d(5:7);
  
    % system states
    q = x(1:4);      % quaternion orientation
    %q = q/norm(q);   % renormalize to stay in S^3, do I really need this?
    w = x(5:7);      % body velocity
    w_quat = [0; w]; % augment velocity for hamilton product

    % drift and actuation matrices
    f = [0.5 * mult(w_quat,q); 
         zeros(3,1)];
    g = [zeros(4,3); eye(3)];

    % compute xi
    e_pos = mult(q_conj(q_d), q);
    e_pos = mult(q, q_conj(q_d));
    xi = Log(e_pos);

    % CLF on lie algebra with (eta) dynamics
    eta = [xi; w];
    V = eta' * P * eta;
    LfV = eta' * (A'*P + P*A) * eta;
    LgV = 2* eta'* P * B; 

    a = LfV + (lam_min_Q/lam_max_P) *V;
    b = LgV';

    % min-norm QP solution (Theorem 25.3)
%     if a > 0
%         u_clf = -0.5 * (LfV / (eta'*P*B*B'*P*eta)) *B'*P*eta;
% %         u_clf = (a*b)/(b'*b)
%     else
%         u_clf = [0;0;0]
%     end

    % LQR solution
    u_clf = -0.5*B'*P*eta;

    % dynamics
    xdot = f + g * (u_clf);

end

% rotating body dynamics
function xdot = quatODE_PD(~, x, x_d, K)
   % desired eq points
   q_d = x_d(1:4);
   w_d = x_d(5:7);
 
   % system states
   q = x(1:4)      % quaternion orientation
   %q = q/norm(q);   % renormalize to stay in S^3, do I really need this?
   w = x(5:7)      % body velocity
   w_quat = [0; w]; % augment velocity for hamilton product

   % drift and actuation matrices
   f = [0.5 * mult(w_quat,q); 
        zeros(3,1)];
   g = [zeros(4,3); eye(3)];

   % error states
   e_pos = mult(q_conj(q_d), q);
   e_vel = w - w_d;
   xi = Log(e_pos);

   % PD control
   K_p = diag([1, 1, 1]);
   K_d = diag([1, 1, 1]);

   u_pd = -K_p*e_pos(2:4) - K_d * e_vel;

   % dynamics
   xdot = f + g * (u_pd);

end

% rotating body dynamics
function xdot = quatODE_inertia(~, x, x_d)
   % desired eq points
   q_d = x_d(1:4);
   w_d = x_d(5:7);

   % sytem paramters
   J = eye(3);
 
   % system states
   q = x(1:4)      % quaternion orientation
   %q = q/norm(q);   % renormalize to stay in S^3, do I really need this?
   w = x(5:7)      % body velocity
   w_quat = [0; w]; % augment velocity for hamilton product

   % drift and actuation matrices
   f = [0.5 * mult(w_quat,q); 
        inv(J)*cross(J*w, w)];
   g = [zeros(4,3); inv(J)];

   % error states
   e_pos = mult(q_conj(q_d), q);
   e_vel = w - w_d;

   % PD control
   K_p = 10* diag([1, 1, 1]);
   K_d = 10* diag([1, 1, 1]);

   u_pd = -K_p*e_pos(2:4) - K_d * e_vel;

   % dynamics
   xdot = f + g * (u_pd);

end

% take quaternion from manifold into the Lie algebra
function elem_lie = Log(q)
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    w = q0;
    v= [q1; q2; q3];

    elem_lie = 2*v*atan2(norm(v),w)/norm(v); % from MicroLie Theory

end

% take the lie group element back into the the manifold
function q = Exp(elem_lie)

    w = elem_lie; % lie algebra elemnt is a pure quaternion

    q = [cos(norm(w)/2);               % real part
         w/norm(w)*sin(norm(w)/2)];    % imaginary parts
end

% take the conjugate of a quaternion
function q_conjugate = q_conj(q)
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    q_conjugate = [q0; -q1; -q2; -q3]; 
end

% quaternion multiplication
function q = mult(q1,q2)

    a = q1(1);
    b = q1(2);
    c = q1(3);
    d = q1(4);

    e = q2(1);
    f = q2(2);
    g = q2(3);
    h = q2(4);

    q(1) = a*e - b*f - c*g - d*h;
    q(2) = a*f + b*e + c*h - d*g;
    q(3) = a*g - b*h + c*e + d*f;
    q(4) = a*h + b*g - c*f + d*e;

    q = q';

end















