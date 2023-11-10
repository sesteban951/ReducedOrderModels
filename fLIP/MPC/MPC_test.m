%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MPC implementaion of simple inverted pendulum system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; clc;

% symbolic time and desired time
syms t t_d     % state = [theta, theta_dot]

% drift and acutation vectors
f = [t_d; -sin(t)];
g = [0; 1];

% Jacobians of f(x) and g(x), function of theta
Df = [diff(f,t) diff(f,t_d)];
Dg = [diff(g,t) diff(g,t_d)];

% turn Df, Dg, f into functions of state, theta, theta_dot
Df = matlabFunction(Df,'Vars',{'t','t_d'});
Dg = matlabFunction(Dg,'Vars',{'t','t_d'});
f = matlabFunction(f, 'Vars',{'t','t_d'});

% define terminal point. In this case, the origin
% also define the inital condition
x_N = [0 0]';
x_0 = [2,3]';

% specify MPC parameters
dt = 0.1;                           % step size
t_span = 10;                          % fwd integration amount
N = 20;                              % number of horizon steps
state_dim = 2;                       % x in R^n 
input_dim = 1;                       % u in R^d
% total_dim = state_dim + input_dim;   % n + d

% weights for cost matrix, each weight
state_scaling = [1 .1];
input_scaling = 1;

% build Qx and Qu matrices for the horizon
state_cost = kron(eye(N),diag(state_scaling));
input_cost = kron(eye(N-1),diag(input_scaling));
Q = blkdiag(state_cost, input_cost);     % [Qx1, ..., QxN, Qu1, ..., QuN-1]
% we do this big matrix to avoid triple indexing / creating  tensor matrix

% size and indeces of vectors
state_size_wo_terminal = state_dim * (N-1); 
% state_size_w_terminal  = state_dim * N;  
input_size = input_dim * (N-1);             % we dont give an inpuit at N

% state_end_idx = state_dim * N;                      % stacked vectors [x;u]
% input_end_idx = state_end_idx + input_dim * (N-1);  % stacked vectors [x;u]

% define optimizer variables
yalmip('clear')                       % good practice top clear this
states = sdpvar(state_dim,N);
inputs = sdpvar(input_dim,N-1);

x0_ = sdpvar(state_dim,1);           % initalize intial cond as opt param

Ad_ = sdpvar(state_dim, state_size_wo_terminal);   % discrete drift
Bd_ = sdpvar(state_dim, input_size);               % discrete actuation
Cd_ = sdpvar(state_dim, N-1);                      % discrete residual

% create containers to hold all linear discrete matrices
Ad = zeros(state_dim, state_size_wo_terminal);
Bd = zeros(state_dim, input_size);
Cd = zeros(state_dim, N-1);

%%%%%%  build constraint list
% input constraints, u in U
u_max = 100;
u_min = -100 ;

% state constraints, x in X
pos_max = 100;
pos_min = -100;
vel_max = 100;
vel_min = -100;

% matrices and vectors to describe polyhedron for state cons., Cx <= d
C = [1, 0; -1, 0; 0, 1; 0, -1];
d = [pos_max; -pos_min; vel_max; -vel_min];

% init constraint list and add the inital constraint to it
init_const = (states(:,1) == x0_) ;
constraints = [init_const];

% for loop to add all constraint
for i = 1:N-1

    % new input constraint
    input_cons = (u_min <= inputs(i)) & (inputs(i) <= u_max);

    % new state constraint, described as polyhedron
    state_cons = (C*states(:,i) <= d);

    % new dynamics constraint
    dyn_cons = states(:,i+1) == Ad_(:,(i-1)*state_dim+1 : i*state_dim)*states(:,i) ... 
               + Bd_(:,i)*inputs(i) ...
               + Cd_(:,i);

    % append to the constraints list
    constraints = [constraints; input_cons; state_cons; dyn_cons];

end

% terminal constraint
term_const = states(:,N) == x_N;
constraints = [constraints; term_const];

% build objective, stack x and u ---> J = [x; u]' * Q * [x;u]
objective = (1/2) * ([states(:); inputs(:)]' * Q *  [states(:); inputs(:)]);

% configure optimizer settings
options = sdpsettings('solver','mosek','verbose',2);

%%%%%% Setup the optimizer (line 185 of MPC_Bez)
P = optimizer(constraints, ...                 % constraints
              objective, ...                   % objective
              options,...                      % options
              {Ad_, Bd_, Cd_, x0_},...         % inpute paramters
              {states,inputs});                % dec. vars.

% init history containers
% x_d = x_N;
T = [0];      %  history of time steps 
X = [x_0'];   %  history of sates at tiem steps
U = [];

% fwd simulate for every time step
for t = 0:dt:t_span
    
    % linear interpolation of current state/input to desired state
    x_bar = [linspace(x_0(1), x_N(1), N); linspace(x_0(2), x_N(2), N)];
    u_bar = [linspace(0,0,N-1)];  % autonamous system

    % Build Ad, Bd, Cd for every integration step
    Ad = [];
    Bd = [];
    Cd = [];

    for k = 1:N-1
        % which x_bar and u_bar
        x_k = x_bar(:,k);
        u_k = u_bar(:,k);
    
        % continous time matrices
        Ac = Df(x_k(1), x_k(2)) + Dg(x_k(1), x_k(2)) * u_k;
        Bc = g;     % g is a constant
        Cc = f(x_k(1), x_k(2)) + g * u_k - Ac * x_k - Bc * u_k;

        % discrete time matrices
        [Ad_k,Bd_k,Cd_k] = css2dss('Exact',dt,Ac,Bc,Cc);

        % append to discrete matrix list
        Ad = [Ad , Ad_k];
        Bd = [Bd , Bd_k];
        Cd = [Cd , Cd_k];

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is where the magic happens -- optimzer plus dynamics solver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % solve MPC problem for x_star and u_star given init state x_0
    % solution will return 
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
    [t,x] = ode45(@(t,x) PendulumODE(t,x,u_ff), [0 dt], x_0);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store time steps and states
    T = [T; t+T(end)]; % add new batch of time steps, add T(end) to
                       % splice and make continous stream of time

    X = [X; x];        % add new batch of states at time steps

    U = [U; ones(length(t), 1) * u_ff]; % add new batch of inputs

    % new inital condition is the last state after MPC u_ff is applied
    x_0 = x(end,:)';

end

% plot states history
plot(T,X)
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('State, $x_i$', Interpreter='latex')
legend('$x_1$','$x_2$',Interpreter='latex')
grid on
drawnow

% plot input history
figure;
plot(T(1:length(T)-1 , :),U)
xlabel('Time, $t$ [sec]', Interpreter='latex')
ylabel('Feed Fwd. Input, $u_{ff}$', Interpreter='latex')
grid on

% function [x_d, u_d] = MPC(t,x)
% %%%% call the optimizer (line 262 of MPC_Bez)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simple function to evaulate dynamics at specific state although the 
% function f at the top kind of does the same thing lol.
function dx = PendulumODE(t,x,u_ff)

    % can construcrt PD control or low level CLF controller here
    % K = [-1 -1];
    % u = K*(x-x_d);

    % only consdier the feedfwd term in the dynamics
    u = u_ff;
    
    % return the dynamics. ode45 usually wants colm. vector
    dx = [x(2); -sin(x(1))] + [0; 1]*u;

end