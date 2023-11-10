%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MPC implementaion of fLIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% system parameters
Ip = 0.5;
If = 0.2;
bp = 0.01;
bf = 0.1;
m1 = 2;
m2 = 1;
L = .5;
g = 9.81;

% store actual system parameters 
Lm = L*0.5;     % CG off set along pole
th_m = pi/4;    % CG off set along pole

params_actual = [Ip;
          If;
          bp;
          bf;
          m1;
          m2;
          L;
          Lm;
          th_m;
          g];

% store wrong system parameters
Lm = L;    % CG off set along pole
th_m = 0;  % CG off set along pole

params_wrong = [Ip;
          If;
          bp;
          bf;
          m1;
          m2;
          L;
          Lm;
          th_m;
          g];

% init function to linearize about x_bar jacobian
[Df_actual, Dg_actual] = jacobians(params_actual);   % <----------------uncertainty
[Df_wrong, Dg_wrong] = jacobians(params_wrong);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define terminal point and the inital condition
x_N = [0; 0; 0; 0];     % flywheel pos doesnt matter
x_0 = [pi; 0; 0; 0];

% specify MPC parameters
dt = 0.05;      % integration step size
t_span = 4;     % total integration time
N = 5;         % number of steps in horizon; T_tot = N * dt

% specify state and input dimension (hardcoded)
state_dim = 4;  % x = [th, th_f, th_dot, th_f_dot]
input_dim = 1;  % u = tau_ext

% state and input size for MPC
state_size = state_dim * (N-1); 
input_size = input_dim * (N-1);
state_size_wo_terminal = state_dim * (N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define optimizer variables
yalmip('clear')                       % good practice to clear this

% super wide vectors as well
states = sdpvar(state_dim,N);      % states, [x1, ..., xN] a (4 X N) matrx
inputs = sdpvar(input_dim,N-1);    % inputs, [u1, ..., uN-1]

x0_ = sdpvar(state_dim,1);           % initalize initial cond as opt param

Ad_ = sdpvar(state_dim, state_size_wo_terminal);   % discrete drift
Bd_ = sdpvar(state_dim, input_size);               % discrete actuation
Cd_ = sdpvar(state_dim, N-1);                      % discrete residual

% containers to hold all linear discrete matrices, super wide matrices
Ad = zeros(state_dim, state_size);
Bd = zeros(state_dim, input_size);
Cd = zeros(state_dim, N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build constraint list
% state constraints, x in X
th_max =  10;           % pole position constraints
th_min = -10;
th_f_max =  1E6;        % flywheel position constraints
th_f_min = -1E6;

th_dot_max =  30;       % pole velocity conastraints
th_dot_min = -30;
th_f_dot_max =  20000;    % flywheel velocity constraints
th_f_dot_min = -20000;

% input constraints, u in U
u_max =  1E2;
u_min = -1E2;

% matrices and vectors to describe polyhedron for state cons., Cx <= d
C = kron(eye(state_dim),[1;-1]);
d = [th_max;     % put all constraints into polyhedron
     -th_min; 
     th_f_max; 
     -th_f_min; 
     th_dot_max; 
     -th_dot_min; 
     th_f_dot_max; 
     -th_f_dot_min];

% init constraint list and add the inital constraint to it
init_const = (states(:,1) == x0_);
constraints = [init_const];

% for loop to add all constraints
for i = 1 : N-1

    % new state constraint, described as polyhedron
    state_cons = (C*states(:,i) <= d);

    % new input constraint
    input_cons = (u_min <= inputs(i)) & (inputs(i) <= u_max);

    % new dynamics constraint
    dyn_cons = states(:,i+1) == Ad_(:,(i-1)*state_dim+1 : i*state_dim)*states(:,i) ... 
               + Bd_(:,i)*inputs(i) ...
               + Cd_(:,i);
    
    % append constraints to the constraints list
    constraints = [constraints; input_cons; state_cons; dyn_cons];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build objective function
% weights for cost matrix, each weight
state_scaling = [1, 0, .1, .1];
input_scaling = .01;

% build Qx and Qu matrices for the horizon in giant Q matrix
state_cost = kron(eye(N),diag(state_scaling));
input_cost = kron(eye(N-1),diag(input_scaling));
Q = blkdiag(state_cost, input_cost);     % [Qx1, ..., QxN, Qu1, ..., QuN-1]

% boolean for whether to use terminal condition or not
term_cond = 1;

% built in terminal constraint into the objective function
if term_cond

    % build objective, stack x, u -> J = [x; u]' * Q * [x;u] + x_N' * V * x_N
    terminal_scaling = [1, 0, 1, 1];
    term_cost = diag(terminal_scaling);
    V = term_cost;

    objective = (1/2)*([states(:); inputs(:)]' * Q * [states(:); inputs(:)])+...
            (states(:,N) - x_N)' * V * (states(:,N) - x_N);

% built in terminal constraint into the constraints
else

    % append terminal constraint -- don't include flywheel pos
    term_const1 = (states(1,N) == x_N(1));
    term_const3 = (states(3,N) == x_N(3));
    term_const4 = (states(4,N) == x_N(4));
    
    constraints = [constraints; term_const1, term_const3, term_const4];

    % build objective, stack x, u -> J = [x; u]' * Q * [x;u] 
    objective = (1/2)*([states(:); inputs(:)]' * Q * [states(:); inputs(:)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the optimizer
% configure optimizer settings
options = sdpsettings('solver','mosek','verbose',0);

% create the optimzer object/function
P = optimizer(constraints, ...                 % constraints
              objective, ...                   % objective
              options,...                      % options
              {Ad_, Bd_, Cd_, x0_},...         % inpute paramters
              {states,inputs});                % dec. vars.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use LQR to compute fixed control policy, u = -K_lqr * x
% linearize about the desired target, x_dot = A_xN x + B_xN u
A_xN = Df_wrong(x_N(1),x_N(3),x_N(4));  % <----------------uncertainty

[~, g_] = dyn_vect(x_N(1),x_N(3),x_N(4), params_wrong); % <----------------uncertainty
B_xN = g_;

% cost matrices
Qx = diag([5 * 1; 0.001; 1; 1]);
Qu = diag(ones(input_dim));

% will use this gain matrix for forward integarting 
K_lqr = lqr(A_xN, B_xN, Qx, Qu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fun stuff is now to run MPC
% init history containers
T = [0];      %  history of time steps 
X = [x_0'];   %  history of states at time steps
U = [];       %  history of inputs at time steps

% contaniers to store MPC trajectories (Noel did this)
X_MPC_Path = zeros (N, state_dim, ceil(t_span/dt));
U_MPC_Path = zeros (N-1, input_dim, ceil(t_span/dt));
T_MPC_Path = zeros(N,1,ceil(t_span/dt));
x_path_ind = 1;

% fwd simulate for every time step
for t = 0:dt:t_span
    
    % linear interpolation of current state/input to desired state
    %{
    x_bar = [linspace(x_0(1), x_N(1), N); 
             linspace(x_0(2), x_N(2), N);
             linspace(x_0(3), x_N(3), N);
             linspace(x_0(4), x_N(4), N)];

    u_bar = [linspace(0,0,N-1)];  % autonamous system dynamics
    %}
    
    % generate trajectory with linearized sys about xN and LQR fbk control
    % fwd integrate until within ball around xN
    %{
    t_lqr = 0:0.01:500;

    opt = odeset('Events', @close_to_xN);
    [t_,~] = ode45(@(t,x) fLIP_ODE_control(t,x,K_lqr, params), ...
                       t_lqr, x_0, opt);

    t_term = t_(end);
    t_lqr = 0:(t_term)/(N-1):t_term;
    [~,x_] = ode45(@(t,x) fLIP_ODE_control(t,x,K_lqr, params), ...
                         t_lqr, x_0, opt);
    x_bar = x_' ;
    u_bar = sum( (-K_lqr.*x_) , 2); %  computing u_k = -K * x_k
    u_bar = u_bar';
    %}

    % fwd simulate for small dt
    t_lqr = 0:dt:dt*(N-1);
    [~,x_] = ode45(@(t,x) fLIP_ODE_control(t,x,K_lqr, params_wrong), ...  % <----------------uncertainty
                         t_lqr, x_0);
    x_bar = x_' ;
    u_bar = sum( (-K_lqr.*x_) , 2); %  computing u_k = -K * x_k
    u_bar = u_bar';

    % Build Ad, Bd, Cd for every integration step
    Ad = [];
    Bd = [];
    Cd = [];

    for k = 1:N-1
        % which x_bar and u_bar
        x_k = x_bar(:,k);
        u_k = u_bar(:,k);

        % compute f(x), g(x) linearizations
        Df = Df_wrong(x_k(1), x_k(3), x_k(4));   % <----------------uncertainty
        Dg = Dg_wrong(x_k(1), x_k(3), x_k(4));   % <----------------uncertainty
        [f_x, g_x] = dyn_vect(x_k(1), x_k(3), x_k(4), params_wrong); % <----------------uncertainty
    
        % continous time linear matrices
        Ac = Df + Dg * u_k;
        Bc = g_x;   
        Cc = f_x + g_x * u_k - Ac * x_k - Bc * u_k;

        % discrete time linear matrices
        [Ad_k,Bd_k,Cd_k] = css2dss('Exact',dt,Ac,Bc,Cc);  %%%%%%%%%%%%%

        % append to discrete matrix list
        Ad = [Ad , Ad_k];
        Bd = [Bd , Bd_k];
        Cd = [Cd , Cd_k];

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is where the magic happens -- optimzer plus dynamics solver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % solve MPC problem for x_star and u_star given init state x_0
    % solution will return optimal trajectory and optimal input
    [sol, diagnostics,d1,d2,d3,d4] = P({Ad,Bd,Cd,x_0});
    x_star = sol{1};    
    u_star = sol{2};

    % take first input and use as feed forward via zero order hold
    u_ff = u_star(1);
    % Note: you can also devise some kind of low level control to track
    % the optimal state, x_star. Can do it here or inside PendulumODE

    % fwd simulate with given x0 and u_ff to see where you end up after the
    % small time step

    % [0,dt] is the span of time in which you want ot integrate
    % x_0 is the intial condition for the solver
    % PendulumODE is the dynamics of your system with inputs 
    [t,x] = ode45(@(t,x) fLIP_ODE(t,x,u_ff,params_actual), [0 dt], x_0);  % <----------------uncertainty

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store time steps and states for actual evolution of system
    t_now = T(end);
    T = [T; t+T(end)]; % add new batch of time steps, add T(end) to
                       % splice and make continous stream of time

    X = [X; x];        % add new batch of states at time steps

    U = [U; ones(length(t), 1) * u_ff]; % add new batch of inputs

    % store MPC optimal inputs and states
    x_star = x_star';

    % store MPC trajectories at each step (Noel did this part)
    X_MPC_Path(:,:,x_path_ind) = x_star;
    U_MPC_Path(:,:,x_path_ind) = u_star;
    T_MPC_Path(:,:,x_path_ind) = (t_now+(0:dt:dt*(N-1)))';
    x_path_ind = x_path_ind+1;

    % new inital condition is the last state after MPC u_ff is applied
    x_0 = x(end,:)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(1)
clf
% plot states history
subplot(2,2,1)
plot(T,X)
title('Actual Dynamics', Interpreter='latex')
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('State, $x_i$', Interpreter='latex')
legend('$\theta$','$\theta_{f}$','$\dot{\theta}$','$\dot{\theta}_{f}$',...
       Interpreter='latex')
grid on

% plot input history
subplot(2,2,2)
plot(T(1:length(T)-1 , :),U);
title('Actual Input', Interpreter='latex')
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('Input, $u$', Interpreter='latex')
grid on

% plot MPC states history
subplot(2,2,3)
hold on;
plot(T,X)
% for i = 1:size(T_MPC_Path,3)
for i = 1:size(T_MPC_Path,3)
    plot(T_MPC_Path(1:end,1,i),X_MPC_Path(1:end,1,i),'k--')
end
title('with MPC Dynamics', Interpreter='latex')
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('State, $x_i$', Interpreter='latex')
legend('$\theta$','$\theta_{f}$','$\dot{\theta}$','$\dot{\theta}_{f}$',...
       Interpreter='latex')
grid on

% plot MPC input history
subplot(2,2,4)
hold on;
plot(T(1:length(T)-1 , :),U);
% for i = 1:size(T_MPC_Path,3)
for i = 1:size(T_MPC_Path,3)
    plot(T_MPC_Path(1:2,:,i),U_MPC_Path(1:2,:,i),'k--')
end
title('with MPC Input', Interpreter='latex')
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('Input, $u$', Interpreter='latex')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxillary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute the drift and actuation vectors of the system
function [f_x, g_x] = dyn_vect(th, th_dot, th_f_dot,params)

    % system parameters
    Ip = params(1);
    If = params(2);
    bp = params(3);
    bf = params(4);
    m1 = params(5);
    m2 = params(6);
    L  = params(7);
    Lm = params(8);   % CG off set along pole
    th_m = params(9);   % CG off set along pole
    g  = params(10);

    % Dynamics, mech. sys.
    D = [(Ip+m2*L*L) If; 
         If          If];
    C = [bp 0; 
         0  bf];
    G = [-m1*g*Lm*sin(th + th_m) - m2*g*L*sin(th);
         0];
    B = [0; 1];
    
    % define q_dot and x_dot
    q_dot = [th_dot; th_f_dot];
    % q_ddot = inv(D) * (-C*q_dot - G + B*u);
    
    % drift and actuation vector for the system
    f_x = [q_dot;     -inv(D)*(C*q_dot + G)];
    g_x = [zeros(2,1); inv(D)*B];

end

% get linearizations of vectors, Df, Dg at point x_k
function [Df, Dg] = jacobians(params)

    syms th th_dot th_f_dot

    % system parameters
    Ip = params(1);
    If = params(2);
    bp = params(3);
    bf = params(4);
    m1 = params(5);
    m2 = params(6);
    L  = params(7);
    Lm = params(8);   % CG off set along pole
    th_m = params(9);   % CG off set along pole
    g  = params(10);
    
    % Dynamics, mech. sys.
    D = [(Ip+m2*L*L) If; 
         If          If];
    C = [bp 0; 
         0  bf];
    G = [-m1*g*Lm*sin(th + th_m) - m2*g*L*sin(th);
         0];
    B = [0; 1];
    
    % define q_dot and x_dot
    q_dot = [th_dot; th_f_dot];
    % q_ddot = inv(D) * (-C*q_dot - G + B*u);
    
    % drift and actuation vector for the system
    f_x = [q_dot;     -inv(D)*(C*q_dot + G)];
    g_x = [zeros(2,1); inv(D)*B];

    % compute jacobians
    D_x2 = zeros(4,1);
    Df = [diff(f_x,th), D_x2, diff(f_x,th_dot), diff(f_x,th_f_dot)];
    Dg = [diff(g_x,th), D_x2, diff(g_x,th_dot), diff(g_x,th_f_dot)];

    % state values to compute
    Df = matlabFunction(Df, 'Vars',{'th','th_dot','th_f_dot'});
    Dg = matlabFunction(Dg, 'Vars',{'th','th_dot','th_f_dot'});


end

% flywheel dynamics differential equation to use in for loop
function x_dot = fLIP_ODE(~,x,u_ff,params)

    % only do feedfoward
    u = u_ff;

    % evaluate drift and actuation ,matricves at specific state
    [f_x, g_x] = dyn_vect(x(1), x(3), x(4),params)   ;

    % return dynamics
    x_dot = f_x + g_x * u;

end

% flywheel dynamics differential equation with feedback control
function x_dot = fLIP_ODE_control(~,x,K,params)
    
    % feedback controller
    u_fb = -K * x;

    % total input, feedabck plus feedforward 
    u =  u_fb;

    % evaluate drift and actuation ,matrices at specific state
    [f_x, g_x] = dyn_vect(x(1), x(3), x(4),params);

    % return dynamics
    x_dot = f_x + g_x * u;

end

% ODE fwd simulate until terminating condition 
% Not the right way of doing this thing
function [val, isterminal, direction] = close_to_xN(~,x)

   % stop simualting once within a ball of xN
   val = (norm(x - zeros(4,1)) <= 0.05);
   isterminal = 1;  % boolean for stopping the iterations
   direction  = 0;  

end



