%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot SRB Traj Opt results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% load data
T = readmatrix('./results/time.csv');
X = readmatrix('./results/states.csv');
U = readmatrix('./results/inputs.csv');

% extract the state components
p = X(:,1:3);   % center of mass pos
q = X(:,4:7);   % orientation quaternion
v = X(:,8:10);  % center of mass vel
w = X(:,11:13); % angular velocity
K = size(T,1);

% plotting mode
mode = 2; % 1: animation, 2: plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animation
if mode == 1

    figure;
    hold on;

    % mesh to use
    mesh = "./meshes/legoman.stl";

    % initialize the pose plot
    dt = T(2) - T(1);
    q0 = quaternion(q(1,1), q(1,2), q(1,3), q(1,4));
    p0 = p(1, :);
    % h = poseplot(q0, p0, ScaleFactor=0.1, MeshFileName=mesh);
    h = poseplot(q0, p0);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    grid on; axis equal; view(3);
    
    % loop through time steps
    k = 1;
    while true

        % get the current pose
        p_  = p(k, :);
        q_  = quaternion(q(k,1), q(k,2), q(k,3), q(k,4));

        % update pose
        set(h, Orientation=q_, Position=p_);

        % set the title
        msg = sprintf('Time: %.2f s', T(k));
        title(msg);

        drawnow;
        
        % pause for a little
        pause(dt);
        
        % update the index
        k = k + 1;
        if k > K
            k = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state and input plots
if mode == 2

    figure;
    hold on;

    % Position
    subplot(4,4,1);
    plot(T, p(:,1), 'r', 'DisplayName', 'p_x');
    xlabel('Time (s)');
    ylabel('Pos X (m)');
    grid on;
    subplot(4,4,5);
    plot(T, p(:,2), 'g', 'DisplayName', 'p_y');
    xlabel('Time (s)');
    ylabel('Pos Y (m)');
    grid on;
    subplot(4,4,9);
    plot(T, p(:,3), 'b', 'DisplayName', 'p_z');
    xlabel('Time (s)');
    ylabel('Pos Z (m)');
    grid on;

    % Velocity
    subplot(4,4,2);
    plot(T, v(:,1), 'r', 'DisplayName', 'v_x');
    xlabel('Time (s)');
    ylabel('Vel X (m/s)');
    grid on;
    subplot(4,4,6);
    plot(T, v(:,2), 'g', 'DisplayName', 'v_y');
    xlabel('Time (s)');
    ylabel('Vel Y (m/s)');
    grid on;
    subplot(4,4,10);
    plot(T, v(:,3), 'b', 'DisplayName', 'v_z');
    xlabel('Time (s)');
    ylabel('Vel Z (m/s)');
    grid on;

    % Orientation (quaternion)
    subplot(4,4,15);
    plot(T, q(:,1), 'm', 'DisplayName', 'q_w');
    xlabel('Time (s)');
    ylabel('Quat W');
    grid on;
    subplot(4,4,3);
    plot(T, q(:,2), 'r', 'DisplayName', 'q_x');
    xlabel('Time (s)');
    ylabel('Quat X');
    grid on;
    subplot(4,4,7);
    plot(T, q(:,3), 'g', 'DisplayName', 'q_y');
    xlabel('Time (s)');
    ylabel('Quat Y');
    grid on;
    subplot(4,4,11);
    plot(T, q(:,4), 'b', 'DisplayName', 'q_z');
    xlabel('Time (s)');
    ylabel('Quat Z');
    grid on;

    % Angular Velocity
    subplot(4,4,4);
    plot(T, w(:,1), 'r', 'DisplayName', 'w_x');
    xlabel('Time (s)');
    ylabel('Ang Vel X');
    grid on;
    subplot(4,4,8);
    plot(T, w(:,2), 'g', 'DisplayName', 'w_y');
    xlabel('Time (s)');
    ylabel('Ang Vel Y');
    grid on;
    subplot(4,4,12);
    plot(T, w(:,3), 'b', 'DisplayName', 'w_z');
    xlabel('Time (s)');
    ylabel('Ang Vel Z');
    grid on;


    figure;
    hold on;

    % input forces
    subplot(3,2,1);
    hold on;
    plot(T(1:end-1), U(:,1), 'r', 'DisplayName', 'F_left_x');
    plot(T(1:end-1), U(:,4), 'r', 'DisplayName', 'F_right_x');
    xlabel('Time (s)');
    ylabel('F_x (N)');
    grid on;
    subplot(3,2,3);
    hold on;
    plot(T(1:end-1), U(:,2), 'g', 'DisplayName', 'F_left_y');
    plot(T(1:end-1), U(:,5), 'g', 'DisplayName', 'F_right_y');
    xlabel('Time (s)');
    ylabel('F_y (N)');
    grid on;
    subplot(3,2,5);
    hold on;
    plot(T(1:end-1), U(:,3), 'b', 'DisplayName', 'F_left_z');
    plot(T(1:end-1), U(:,6), 'b', 'DisplayName', 'F_right_z');
    xlabel('Time (s)');
    ylabel('F_z (N)');
    grid on;

    % input moments
    subplot(3,2,2);
    hold on;
    plot(T(1:end-1), U(:,7), 'r', 'DisplayName', 'M_left_x');
    plot(T(1:end-1), U(:,10), 'r', 'DisplayName', 'M_right_x');
    xlabel('Time (s)');
    ylabel('M_x (N.m)');
    grid on;
    subplot(3,2,4);
    hold on;
    plot(T(1:end-1), U(:,8), 'g', 'DisplayName', 'M_left_y');
    plot(T(1:end-1), U(:,11), 'g', 'DisplayName', 'M_right_y');
    xlabel('Time (s)');
    ylabel('M_y (N.m)');
    grid on;
    subplot(3,2,6);
    hold on;
    plot(T(1:end-1), U(:,9), 'b', 'DisplayName', 'M_left_z');
    plot(T(1:end-1), U(:,12), 'b', 'DisplayName', 'M_right_z');
    xlabel('Time (s)');
    ylabel('M_z (N.m)');
    grid on;

end