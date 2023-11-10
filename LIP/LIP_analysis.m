%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To gain intuition behind LIP model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIP model parameters
g = 9.81;

dz = 0.05;
z0 = 0.1:dz:2;

lam = sqrt(g./z0);

% plot different lambda parameters
figure(1);
subplot(1,2,1)
plot(z0,lam,'b')
xlabel("$z_0$, Constant Height",'Interpreter','latex')
ylabel("$\lambda = \sqrt{g/z_0}$, ",'Interpreter','latex')
grid on; 

% find the poles of the drift matrix;
subplot(1,2,2)
xline(0); yline(0);
hold on; grid on; 
xlabel("$z_0$, Constant Height",'Interpreter','latex')
ylabel("$\lambda$, Re-Axis",'Interpreter','latex')
for i = 1:length(lam)
    
    % build drift matrix and compute poles, its just lambda on the real axis
    A = [0,     1;
     lam(i)^2, 0];
    poles = eigs(A);
    plot([z0(i),z0(i)],poles','xr')

end

% consensus is that the theory of having a taller pendulum is confirmed. 