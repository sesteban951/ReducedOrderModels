%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fLIP CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = fLIP_control(x,params)

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

% get state space
th = x(1);
th_f = x(2);              % position of the flywheel doesnt affect dynamics
th_dot = x(3);
th_f_dot = x(4);

% basic PD control
th_ref = 0;
p = 60;
d = 6.2;
u = p*(th-th_ref) + d*(th_dot) ;  % PD-control

end