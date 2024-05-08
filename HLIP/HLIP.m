%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIP, Point Mass Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all;

% LIP parameters
params.g = 9.81;                % gravity
params.z0 = 0.5;                  % constant height
lam = sqrt(params.g/params.z0); % lambda

% gait parameters
params.T = 0.4;                  % single support phase length
params.z_apex = params.z0*0.25;  % peak foot height
params.A = [0     1;
            lam^2 0];            % ssp drift matrix

% x = [pos, vel], i.e., [(CoM - foot_stance), (CoM vel)]
x0 = [0.1;-0.1];

% simulation parameters
dt = 0.02;
tspan = 0:dt:params.T;
n_steps = 10;

% containers to hold data
X = zeros(n_steps*length(tspan),length(x0)); % state
B = zeros(n_steps*length(tspan),length(x0)); % Bezier
U = zeros(n_steps*length(tspan),1);          % foot placement inputs
T = zeros(1,n_steps*length(tspan));          % time

% simulate a step
x_current = x0;
t_current = 0;
for j = 1:n_steps
    % simulate one swing phase
    for i = 1:length(tspan)
 
        % get state and input at each time
        x = LIP_ssp(tspan(i),x_current,params);
        u = compute_pf(x,params);

        % compute swing foot trajectory
        Px = [0, 0, 0, u/2, u, u , u];
        Pz = [0, 0, 0, (8/3)*params.z_apex, 0, 0, 0]; % 16/5 based on analysis
        b = [bezier(tspan(i)/params.T,Px); 
             bezier(tspan(i)/params.T,Pz)];
        
        % store into containers
        X(i+((j-1)*length(tspan)),:) = x';
        B(i+((j-1)*length(tspan)),:) = b';
        U(i+((j-1)*length(tspan))) = u;
    end

    % update time array
    T((j-1)*length(tspan)+1 : (j)*length(tspan)) = tspan + t_current;

    % current state and time
    x_current = LIP_reset(x,u);
    t_current = tspan(end)*j;

end

% animate the gait
time_scale = 0.75;
animate(X,B,U,T,params,time_scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILLARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solution to linear inverted pendulum in single support
function x_t = LIP_ssp(t,x0,params)
    % compute the closed form solution, x(t) = x0*e^(At)
    x_t = expm(params.A*t) * x0;
end

% function reset dynamics, switch feet when foot strikes the ground
function x_t = LIP_reset(x,u)
    % unpack reset map data
    p_stance = u;
    p_CoM = x(1);
    v_CoM = x(2);

    % reassign state, this is the reset map
    p = p_CoM - p_stance;
    v = v_CoM;
    x_t = [p;
           v];
end

% solution compute foot location placement, simple Raibert foot placement
function p_f = compute_pf(x,params)
    % unpack state info
    pos = x(1);
    vel = x(2);

    % raibert gains
    a = params.alpha;
    b = params.beta;

    p_f = a*vel + b*pos;
end

% compute current foot location bezier
function b = bezier(t,P)
    % compute bezier curve
    n = length(P);
    b = 0;
    for i = 1:n
        b = b + coeff(t, i - 1, n - 1) * P(i);
    end
end

% get coefficients of bezier polynomials
function c = coeff(t, i, n)
    c = nchoosek(n, i) * t^i * (1 - t)^(n - i);
end

% animate the LIP
function animate(X,B,U,tspan,params,time_scale)
    % figure settings
    figure;
    xline(0); yline(0);
    ylabel("Vertical Pos., $z$",'Interpreter','latex','FontSize',15)
    xlabel("Horizontal Pos., $(x)$",'Interpreter','latex','FontSize',15)
    xlim([min(U(:,1))*1.2, max(U(:,1))*1.2]);
    ylim([0, params.z0*1.2]);
    axis equal; grid on; hold on;
%     set(gcf, 'Position', get(0, 'Screensize'));
    
    % animate the plot
    tic
    t_now = tspan(1);
    i = 0;
    tspan = (1/time_scale)*tspan;
    while t_now < tspan(end)
        
        % break if went passed max tspan
        if i==length(tspan)
            break
        end

        % find current time
        t_now = toc;
        i = find(tspan>t_now,1);
        
        % draw legs
        sl = plot([0,X(i,1)],[0,params.z0],'b');      % stance leg
        sw = plot([X(i,1),B(i,1)],[params.z0,B(i,2)],'--','Color',[.1 .1 .1]); % swing leg

        % draw points
        CoM = plot(X(i,1),params.z0,'.k','MarkerSize',50); % CoM
        sf =  plot(0,0,'.b','MarkerSize',20);              % stance foot
        b  =  plot(B(i,1),B(i,2),'.k','MarkerSize',20);    % swing  foot
        fpt = plot(U(i),0,'xr','MarkerSize',20);           % foot placement target
        
        txt = sprintf("Time: %.3f \n $p =$  %.3f \n $v =$ %.3f", ...
                      tspan(i)*time_scale,X(i,1),X(i,2));
        title(txt,'Interpreter','latex','FontSize',15)
        
        drawnow

        % remove stuff from from plot
        cleanPlot = [sl,sw,CoM, sf, fpt,b];
        
        % do not clean plot if last time point
        if i~=length(tspan)
            delete(cleanPlot);
        end
    end
end


