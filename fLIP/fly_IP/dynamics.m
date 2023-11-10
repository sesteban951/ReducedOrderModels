% place to compute all your dynamics

syms Ip If m1 m2 l L g th psi th_dot psi_dot bp bf alph u

% machaincal sysrtem dynamics
D = [Ip + m1* l^2 + If + m2*L^2 , If;
     If                         , If]

D_inv = inv(D)

C = [bp, 0;
     0 , bf]

G = [-m1*g*l*sin(th + alph) - m2*g*sin(th);
      0]

B = [0;
     1]

% configuration space
q_dot = [th_dot; psi_dot];

% full control affine form system dynamics
f = [q_dot;
     D_inv * (-C*q_dot -G)]

g = [zeros(2,1); D_inv*B]

xdot = f + g * u