clear all; clc;

L = [.1503; .320; .129];
B = [0 0.5];

h = Harpy(B,L);

q1 = pi/8; % hip pitch
q2 = pi/3; % 
q = [q1;q2]

x = [0.1;0]

q0 = [pi/7, pi/4]';

q = h.inv_kin(x(:,end),q0)
x = h.fwd_kin(q)

plot(x(1,:),x(2,:)); axis equal; grid on;
xline(0); yline(0)