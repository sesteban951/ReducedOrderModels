%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALIP simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; 

% system parameters
params.m = 1.0;   % mass
params.g = 9.81;  % gravity
params.z0 = 0.5;  % center of mass height

% simualtion params
params.dt = 0.01; % time step

% initial state
x0 = [0.1;   % p, position of the COM relative to stance
      0.1];  % L, angular momentum about stance

% simulation time
tmax = 1.0;
tspan = 0:params.dt:tmax;

% forward simulation
[t, x] = ode45(@(t,x) alip_dynamics(t, x, params), tspan, x0);

% simulate the dynamics
plot(x(:,1),x(:,2),'r-','LineWidth',2); hold on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute a solution to the ALIP dynamics
function xdot = alip_dynamics(t, x, params)

    % ALIP continuous dynamics matrix
    m = params.m;
    g = params.g;
    z0 = params.z0;

    % dynamics matrices
    A = [0,   1/(m*z0);
         m*g, 0];
    B = [0; 1];

    % compute the ankle torque
    u = alip_control(t, x, params);

    % compute the dynamics
    xdot = A*x + B*u;
end

% compute the control action
function u = alip_control(t, x, params)
    
    % current state
    p = x(1); % position of the COM relative to stance
    L = x(2); % angular momentum about stance
    v = L/(params.m *params.z0); % velocity of the COM relative to stance
    
    % ALIP control law
    kp = 100.0; % position gain
    kv = 1.0;  % velocity gain
    u = kp * (0 - p) + kv * (0 - v);
end