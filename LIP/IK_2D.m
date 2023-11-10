%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harpy Inverse Kinematics (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% desired configuration
q1 = pi/8; % hip pitch
q2 = pi/3; % 

% link lengths in [m]
params.L1 = 150.261/1000;
params.L2 = 320/1000;
params.L3 = 129/1000;

% plot answer
p = fwd_kin(q1,q2,params)
subplot(1,2,1)
plot_fwd(q1,q2,params)

% test ik
q = inv_kin(p(1),p(2),params,[pi/4;pi/4]);
subplot(1,2,2)
plot_fwd(q(1),(q2),params)

% compute IK
q_ans = [q1;q2]
q

% compute forward kinematics
function p = fwd_kin(q1,q2,params)
    L1 = params.L1;
    L2 = params.L2;
    L3 = params.L3;

    x =  (L1+L3)*cos(q1) + L2*cos(q1+q2) + L3*cos(q1);
    z = -(L1+L3)*sin(q1) - L2*sin(q1+q2) - L3*sin(q1);

    p = [x;z];
end

% compute inverse kinematics, Newton raphson method
function q = inv_kin(x,z,params,qk)
    
    % IK settings
    err_tol = 0.001;
    max_iter = 500;

    % iterate until convergence
    err = Inf;
    counter = 0;
    while (err > err_tol)
        % compute error so far
        e = [x;z] - fwd_kin(qk(1),qk(2),params);

        % jacobian of old guess
        J = jacobian(qk(1),qk(2),params)

        % iterate for a new guess
        qk = qk + inv(J) * e;

        err = norm(e)
        counter = counter+1
        if counter > max_iter
            break
        end
    end
    
    q =qk;

end

% compute jacobian
function J = jacobian(q1,q2,params)
    L1 = params.L1;
    L2 = params.L2;
    L3 = params.L3;

    J(1,1) = -(L1+L3)*sin(q1) - L2*sin(q1+q2);
    J(1,2) = -L2*sin(q1+q2);
    J(2,1) = -(L1+L3)*cos(q1) - L2*cos(q1+q2);
    J(2,2) = -L2*cos(q1+q2);
end

% plot forward kin
function plot_fwd(q1,q2,params)

    L1 = params.L1;
    L2 = params.L2;
    L3 = params.L3;

    p_hip = [0; 0];
    p_knee = [L1*cos(q1); -L1*sin(q1)];
    p_ankle = p_knee + [L2*cos(q1 + q2); -L2*sin(q1+q2)];
    p_foot = p_ankle + [L3*cos(q1); -L3*sin(q1)];
    P = [p_hip, p_knee, p_ankle,p_foot];

    hold on; grid on; axis equal;
    xline(0); yline(0)
    plot(P(1,:),P(2,:),'k','LineWidth',2)
    plot(P(1,:),P(2,:),'.b','MarkerSize',20)
    plot(P(1,end),P(2,end),'.r','MarkerSize',25)
    xlabel("x, [m]"); ylabel("z, [m]")

end