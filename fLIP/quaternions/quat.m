%  messing around with quaternions and their Lie Algebras

clear all; clc

% generate random quaternion in S^3
q = rand_quat()

% initial condition
x_old = [1, 0, 0]'

% desired rotations
roll = 0 
pitch = pi
yaw = pi/2

fprintf("-------------------------  Rot  ------------------------------\n")
R_0 = rot_matrix(roll,pitch,yaw) % wrt to world frame (extrinsic)

x_rot = R_0*x_old

fprintf("-------------------------  Quat  -----------------------------\n")
q_0 = rot2quat(R_0)

x_old_quat = to_quat(x_old);
x_quat = mult(q_0, mult(x_old_quat,q_inv(q_0))); % quaternions are associative
x_quat = x_quat(2:4)

fprintf("-------------------------  Q1  ------------------------------\n")
% Q1 log(pq) = log(p) + log(q)? No.
q1 = rand_quat()
q2 = rand_quat()

ans0 = mult(q_conj(q1),q2)
ans1 = Log(ans0)
ans2 = Log(q_conj(q1)) + Log(q2)

% it's not true at all. look at Baker-Campbell-Hausdorff theorem
% see Quaternion Algebra and Calculus by David eberly

%%

q0 = [0, 0, 0, 1];
q0 = q0/norm(q0);
q = quaternion(q0(1), q0(2), q0(3), q0(4))

figure(3);
poseplot(q)

q0 = -[.707, 0, 0, .707];
q0 = q0/norm(q0);
q = quaternion(q0(1), q0(2), q0(3), q0(4))

figure(3);
poseplot(q)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate random quaternion
function q = rand_quat()

    rand_R4 = rand(4,1);
    q = rand_R4 - mean(rand_R4);
    q = (q) / norm(q);

end

% decompose the quaternion into vector to rotate about and the angle of rot
function [alpha, u] = q_decomp(q)
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    % q = cos(alph/2) + (ux*sin(alph/2) + (uy*sin(alph/2) + (uz*sin(alph/2)
    alpha =  2 * acos(q0);
    
    ux = q1/sin(alpha/2);
    uy = q2/sin(alpha/2);
    uz = q3/sin(alpha/2);

    u = [ux; uy; uz];

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

    % wth, dont need the second case for some reason. Maybe because the
    % atan2 handles negative values?
%     % consider the double wraping of lie lagebra with SO(3)
%     if (q0 >= 0)
%         q0 = q(1);
%         q1 = q(2);
%         q2 = q(3);
%         q3 = q(4);
% 
%         w = q0;
%         v = [q1; q2; q3];
% 
%         elem_lie = 2*v*atan2(norm(v),w)/norm(v); % from MicroLie Theory
%     else
%         q0 = -q(1);
%         q1 = -q(2);
%         q2 = -q(3);
%         q3 = -q(4);
% 
%         w = q0;
%         v = [q1; q2; q3];
% 
%         % lie algebra elemnt is a pure quaternion
%         elem_lie = 2*v*atan2(norm(v),w)/norm(v); % from MicroLie Theory
%     end
end

% take the lie group element back into the the manifold
function q = Exp(elem_lie)

    w = elem_lie; % lie algebra elemnt is a pure quaternion

    q = [cos(norm(w)/2);               % real part
         w/norm(w)*sin(norm(w)/2)];    % imaginary parts
end

% convert quaternion to roation matrix
function  rot = quat2Rot(q)

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    r11 = q0^2 + q1^2 - q2^2 - q3^2;
    r12 = 2 * (q1*q2 - q1*q3);
    r13 = 2 * (q1*q3 + q0*q2);
    r21 = 2 * (q1*q2 + q0*q3);
    r22 = q0^2 - q1^2 + q2^2 - q3^2;
    r23 = 2 * (q2*q3 - q0*q1);
    r31 = 2 * (q1*q3 - q0*q2);
    r32 = 2 * (q2*q3 + q0*q1);
    r33 = q0^2 - q1^2 - q2^2 - q3^2;

    rot = [r11, r12, r13;
           r21, r22, r23;
           r31, r32, r33];
end

function q = rot2quat(R)

    r11 = R(1,1);
    r12 = R(1,2);
    r13 = R(1,3);
    r21 = R(2,1);
    r22 = R(2,2);
    r23 = R(2,3);
    r31 = R(3,1);
    r32 = R(3,2);
    r33 = R(3,3);

    qw = sqrt(0.25*(1 + r11 + r22 + r33));
    qx = sqrt(0.25*(1 + r11 - r22 - r33));
    qy = sqrt(0.25*(1 - r11 + r22 - r33));
    qz = sqrt(0.25*(1 - r11 - r22 + r33));

    q = [qw; qx; qy; qz];

end

% take the conjugate of a quaternion
function q_conjugate = q_conj(q)
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    q_conjugate = [q0; -q1; -q2; -q3]; 
end

% take the inverse of a quaternion
function q_inverse = q_inv(q)

    q_inverse = q_conj(q) / (norm(q)^2);
    
end
    
% convert R^3 vector to quat
function q = to_quat(v)

    q = [0 ; v];

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

% roll pitch yaw matrix in that order
function R = rot_matrix(r,p,y)

    R_r = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
    R_p = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
    R_y = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
    
    R = R_y * R_p * R_r;

end
