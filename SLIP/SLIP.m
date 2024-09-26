%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% SLIP params
params.k = 150;   % spring constant
params.br = 0.1;  % damping along spring
params.bt = 0.1;  % damping along pivot point
params.m = 22;    % CoM mass
params.g = 9.81;  % gravity
params.l0 = 0.5;  % spring free length

% make dynamics functions
[f_ground, f_flight] = make_func(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sim params
dt = 0.01;
t_max = 10.0;
t_inf = 10E3;
tspan = [0, t_inf];  % to allow for switching before timeout

% set the switching conditions
opt_g2f = odeset('Events', @(t,x)guard_g2f(t,x,params));
opt_f2g = odeset('Events', @(t,x)guard_f2g(t,x,params));

% set intial condition and which domain
in_flight = true;
if in_flight == true
    % flight: x = [x, z, x_dot, z_dot]
    x0 = [0.0;
          2.0;
          0.1;
          0.0];
    [t,x] = ode45(@(t,x)dynamics_f(t,x,f_flight),tspan,x0,opt_f2g);
elseif in_flight == false
    % ground: x = [r, theta, r_dot, theta_dot]
    x0 = [1.0;
          pi/2;
          0.0;
          0.1];
    [t,x] = ode45(@(t,x)dynamics_g(t,x,f_ground),tspan,x0,opt_g2f);
end

T = [t];
X = [x];
t_current = t(end);
x_current = x(end,:);
while t_current <= t_max

    in_flight = not(in_flight)
    if in_flight == true
        [t,x] = ode45(@(t,x)dynamics_f(t,x,f_flight),[0:dt:t_inf],x_current,opt_f2g);
    elseif in_flight == false
        [t,x] = ode45(@(t,x)dynamics_g(t,x,f_ground),[0:dt:t_inf],x_current,opt_g2f);
    end
    fprintf("-------------------------------")
    T = [T;t+t_current];
    X = [X;to_cartesian(x)];
    t_current = t(end);
    x_current = x(end,:);

end

plot(X(:,1),X(:,2)); hold on;
plot(X(1,1),X(1,2),'.g');
plot(X(end,1),X(end,2),'.r');
xline(0); yline(0);

% % plot CoM location
% cart = to_cartesian(x);

% figure(1); axis equal;
% drawSLIP(t,cart,params);

% figure(2);
% set(gcf,'position',[0,0,700,700])
% subplot(1,2,1);
% plot(cart(:,1),cart(:,2),'b'); hold on;
% plot(cart(1,1),cart(1,2), '.g');
% plot(cart(end,1),cart(end,2), '.r');
% xline(0); yline(0);
% title("CoM Position, $[x ,z]$",'Interpreter','latex')
% xlabel("$x$",'Interpreter','latex')
% ylabel("$z$",'Interpreter','latex')
% grid on; axis equal;

% subplot(1,2,2);
% plot(cart(:,3),cart(:,4),'b'); hold on;
% plot(cart(1,3),cart(1,4), '.g');
% plot(cart(end,3),cart(end,4), '.r');
% xline(0); yline(0);
% title("CoM Position, $[x ,z]$",'Interpreter','latex')
% xlabel("$x$",'Interpreter','latex')
% ylabel("$z$",'Interpreter','latex')
% grid on; axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% forward propagate the dynamcis
function xdot = dynamics_f(t,x,f_flight)

    xdot = f_flight(x(1),x(2),x(3),x(4));

end
function xdot = dynamics_g(t,x,f_ground)
    
    xdot = f_ground(x(1),x(2),x(3),x(4));

end

% Switching surface, ground --> flight
function [value, isterminal, direction] = guard_g2f(t,x,params,tmax)
    
    l0 = params.l0;
    r = x(1);

    value = (r >= l0); % g --> f: r>= l0
    isterminal = 1;
    direction = 0;
end

% Switching surface, flight --> ground
function [value, isterminal, direction] = guard_f2g(t,x,params)
    
    l0 = params.l0;
    r = x(1);
    th = x(2);
    z = r*cos(th);

    value = (z <= l0*cos(th));   % f --> g: z <= l0*cos(th)
    isterminal = 1;
    direction = 0;
end

% polar to cartersian coordinates
function cart = to_cartesian(x)
    % unpack states
    r = x(:,1);
    th = x(:,2);
    rdot = x(:,3);
    thdot = x(:,4);
    
    % convert to cartesian coordinates
    cart = zeros(length(r),4);
    for i = 1:length(r)
        pos_x = -r(i)*sin(th(i));
        pos_z =  r(i)*cos(th(i));
        vel_x =  -rdot(i)*sin(th(i)) - r(i)*cos(th(i))*thdot(i);
        vel_z =   rdot(i)*cos(th(i)) - r(i)*sin(th(i))*thdot(i);
        cart(i,:) = [pos_x, pos_z, vel_x, vel_z];
    end
end

% functions to evaluate SLIP dynamics
function [f_ground, f_flight] = make_func(params)
    
    % states
    syms r theta r_dot theta_dot % ground
    syms x z x_dot z_dot         % flight

    % params
    k = params.k;
    br = params.br;
    bt = params.bt;
    m = params.m;
    g = params.g;
    l0 = params.l0;

    % ground dynamics, x = [r, theta, r_dot, theta_dot]
    f_g = [r_dot;
           theta_dot;
           r*theta_dot^2 - g*cos(theta) + (k/m)*(l0-r) - (br/m);
           -(2/r)*r_dot*theta_dot + (g/r)*sin(theta) - bt/(m*r^2)*theta_dot];
    % flight dynamics, x = [x, z, x_dot, z_dot]
    f_f = [x_dot;
           z_dot;
           0;
           -g];
    
    % make matlab functions to eval the dynamics
    f_ground = matlabFunction(f_g,'Vars',{'r','theta','r_dot','theta_dot'});
    f_flight = matlabFunction(f_f,'Vars',{'x','z','x_dot','z_dot'});
end

