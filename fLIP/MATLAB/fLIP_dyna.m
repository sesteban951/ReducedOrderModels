%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fLIP DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_dot = fLIP_dyna(t,x,params)

% system parameters
Ip = params(1);
If = params(2);
bp = params(3);
bf = params(4);
m1 = params(5);
m2 = params(6);
L = params(7);
Lm = params(8);
g = params(9);
th_m = params(10);

% state space variables
th = x(1);
th_f = x(2);              % position of the flywheel doesnt affect dynamics
th_dot = x(3);
th_f_dot = x(4);

% Dynamics, mech. sys.
D = [(Ip+m2*L*L) If; 
     If          If];
C = [bp 0; 
     0  bf];
G = [-m1*g*Lm*sin(th + th_m) - m2*g*L*sin(th);
     0];
B = [ 0;
      1];

% define control policy
tau_ext = fLIP_control(x,params);
u = tau_ext;

% define q_dot and x_dot
q_dot = [th_dot; th_f_dot];
q_ddot = inv(D) * (-C*q_dot - G + B*u);

f_x = [q_dot;     -inv(D)*(C*q_dot + G)];
g_x = [zeros(2,1); inv(D)*B];

x_dot(1) = q_dot(1);
x_dot(2) = q_dot(2);
x_dot(3) = q_ddot(1);
x_dot(4) = q_ddot(2);

% return xdot
x_dot = x_dot';
