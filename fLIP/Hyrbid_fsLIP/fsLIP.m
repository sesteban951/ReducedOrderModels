%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid sfLIP (2d hoppinf robot with flywheel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; clf; clear figures;


%model parameters
params.g = 9.81;
params.m = 1;
params.Ip = 1;
params.If = 1;
params.k = 1;

% inital and desired states
x0 = [0;    
      0;    
      0;
      0;
      0;
      0;    
      0;    
      0;
      0;
      0];

% simulate dynamics
hz = 45;
t_resolution = 1/hz;
tmax = 10;
t_span = 0 : t_resolution : tmax;

[t,x] = ode45(@(t,x)dynamics(t,x,params),t_span,x0);

% plot some stuff



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dynamics
function xdot = dynamics(t,x,params)

    % states
    y = x(1);
    z = x(2);
    thp = x(3);
    thf = x(4);
    l = x(5);
    ydot = x(6);
    zdot = x(7);
    thpdot = x(8);
    thfdot = x(9);
    ldot = x(10);

    % mechancial matrices
    D = [params.m 0 0         0         0;
         0 params.m 0         0         0;
         0 0 params.Ip 0         0;
         0 0 0         params.If 0;
         0 0 0         0         params.m];
    H = [0;
         -params.m * params.g;
         params.m * params.g * l * sin(thp);
         0;
         -params.k * l];
    B = [0, 0;
         0, 0;
         0, 0;
         1, 0;
         0, 1];

    % configuration vars
    qdot = [ydot;zdot;thpdot;thfdot;ldot];

    % get controller
    u = controller(t,x,params);

    % dynamics
    xdot = zeros(length(x),1);
    xdot(1:length(qdot)) = qdot;
    xdot(length(qdot)+1:length(x)) = inv(D)*(-H+B*u);

end

% controller
function u = controller(t,x,params)

    u = [0;0];

end
