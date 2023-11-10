% 
clear all; clc; close all;

% check that Log and Exp are reliable
% for i = 1:10000
%     fprintf('------------------------------------------------')
%     i
%     r = randrot();
%     [qw, qx, qy, qz] = parts(r);
%     q = [qw, qx, qy, qz]
%     s = Log(q)'
%     e = Exp(s')'
% 
%     err = norm(e-q)
% 
%     if err>0.00001
%         fprintf("found error")
%     end
% % end

% interpolate from random quat to upside down
% r = randrot();
% [qw, qx, qy, qz] = parts(r);
% q0 = [qw, qx, qy, qz];
% s0 = Log(q0);
% s0=[0,0,0]';
% 
% q_mid = [0,1,0,0];
% s_mid = Log(q_mid);
% 
% t = 0:0.005:1;
% 
% st1 = zeros(3,length(t));
% qt1 = zeros(4,length(t));
% for i = 1:length(t)
%     st_ = (1-t(i)) * s0 + (t(i)) *s_mid;
%     st1(:,i) = st_;
% 
%     qt_ = Exp(st_);
%     qt1(:,i) = qt_;
% end
% 
% %interpolate from upsode down to upright
% qf = [-0.999,0,0,0];
% sf = Log(qf);
% sf = [0,0,0]';
% 
% st2 = zeros(3,length(t));
% qt2 = zeros(4,length(t));
% for i = 1:length(t)
%     st_ = (1-t(i)) * s_mid + (t(i)) *sf;
%     st2(:,i) = st_;
% 
%     qt_ = Exp(st_);
%     qt2(:,i) = qt_;
% end

% qt = [qt1, qt2];
% [rows, cols] = size(qt);
% for i = 1:cols
%     q = quaternion(qt(1,i), qt(2,i), qt(3,i), qt(4,i));
%     poseplot(q)
%     drawnow
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using SLERP (spherical linear interpolation)

% t = 0:0.01:1;
% 
% r = randrot();
% [qw, qx, qy, qz] = parts(r);
% q1 = [qw, qx, qy, qz];
% 
% % doing it with two separate interpolations is bad.
% q2 = [0,1,0,0];
% q2 = q2/norm(q2);
% 
% q = [];
% for i = 1:length(t)
%     q_ = slerp(q1,q2,t(i));
%     q = [q q_'];
% end
% 
% q3 = [1,0,0,0];
% q3 = q3/norm(q3);
% for i = 1:length(t)
%     q_ = slerp(q2,q3,t(i));
%     q = [q q_'];
% end

% for i =1:length(q)
%     q_ = q(:,i);
%     poseplot( quaternion(q_(1),q_(2),q_(3),q_(4)) );
%     drawnow
% end

% do it with one interpolation but take the long geodesic
% % instead of the short one
% q3 = [0.707,0.707,0,0];
% q = [];
% for j = 1:length(t)
%     q_ = slerp(q1,q3,t(j));
%     q = [q, q_'];
% end
% 
% for i =1:length(q)
%     q_ = q(:,i);
%     poseplot( quaternion(q_(1),q_(2),q_(3),q_(4)), MeshFileName="torso_reduced.stl", ScaleFactor=0.1 );
%     drawnow
% end

% %% using rodriguez formula to rotate about axis
% 
% th = 0:0.05:2*pi;
% 
% % r = randrot();
% % [qw, qx, qy, qz] = parts(r);
% % q1 = [qw, qx, qy, qz];
% q1 = [1,0,0,0];
% R1 = quat2rotm(q1);
% 
% p0 = [0,0,0];
% pf = [1,1,0];
% z_dir = -[0,0,1];
% 
% u = cross(z_dir,(pf-p0));
% u = u/norm(u);
% S = toSkew(u);
% 
% subplot(1,2,1)
% poseplot(R1, MeshFileName="batman2_repaired.stl", ScaleFactor=0.001);
% for i = 1:length(th)
%     R2 = computeR(th(i),S);
%     p = ((1-th(i))/th(end))*p0 + (th(i)/th(end))*pf;
%     subplot(1,2,2)
%     q = rotm2quat(R2*R1);
%     q = quaternion(q(1),q(2), q(3), q(4));
%     poseplot(q, p, MeshFileName="batman2_repaired.stl", ScaleFactor=0.001);
%     xline(p0(1)); xline(pf(1)); 
%     xline(p0(2)); yline(pf(2)); 
%     drawnow
% end

%% using quaternions

% get random quaternion
r = 2*rand(2,1)-1;
q0 = [r(1),0,0,r(2)];
% q0 = [0,0,0,1];
q0 = q0/norm(q0);
q = quaternion(q0(1),q0(2), q0(3), q0(4));

% positions
p0 = [0,0,0];
pf = [1,1,0];
zhat = [0,0,-1];
v = cross(zhat,pf-p0);
v = v/norm(v);

subplot(2,1,1)
poseplot(q, MeshFileName="batman2_repaired.stl", ScaleFactor=0.001);
% angle = 2*quat_dist(q0,[1,0,0,0])
angle = 2*pi;
th = 0:0.1:angle;
for i = 1:length(th)
    
    q_des(1) = cos(th(i)/2);
    q_des(2) = v(1)*sin(th(i)/2);
    q_des(3) = v(2)*sin(th(i)/2);
    q_des(4) = v(3)*sin(th(i)/2);
    q_interp = mult_quat(q_des,q0);
    q = quaternion(q_interp(1),q_interp(2), q_interp(3), q_interp(4));

    p = ((1-th(i))/th(end))*p0 + (th(i)/th(end))*pf;
    
    subplot(2,1,2)
    poseplot(q, p, MeshFileName="batman2_repaired.stl", ScaleFactor=0.001);
    xline(p0(1)); xline(pf(1)); 
    xline(p0(2)); yline(pf(2)); 
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% geodesic distance
function dist = quat_dist(q1,q2)

    q1_ = [q1(1),q1(2),q1(3),q1(4)];
    q2_ = [q2(1),q2(2),q2(3),q2(4)];

    % opposite hemispehre
    if (q1_*q2_' >= 0)
        q2_ = -q2_;
    end

%     dist = acos(2*(q1_*q2_')^2-1);
    dist = acos(q1_*q2_');
end

% spherical linear interpolation
% https://stackoverflow.com/questions/62943083/interpolate-between-two-quaternions-the-long-way
function q = slerp(q1, q2, t) 
    
    angle = acos(q1 * q2');
    denom = sin(angle);
    
    % check if denom is zero
    q = (q1*sin((1-t)*angle)+q2*sin(t*angle))/denom;
end

% compute rodriguez
function R = computeR(th,Omega)
    R = eye(3) + sin(th) * Omega + (1- cos(th))*(Omega^2);
end

% axis of rotation
function Omega = toSkew(u)
    
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);

    Omega = [0, -u3, u2;
             u3, 0, -u1;
             -u2, u1, 0];

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

% quaternion multiplication
function q = mult_quat(q1,q2)

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