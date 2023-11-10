%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIM fLIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; close all; clc;

% system parameters
Ip = 0.5;
If = 0.5;
bp = 0.5;
bf = 0.05;
m1 = 1;
m2 = 1;
L = 1;%
Lm = L*1;    % CG off set along pole
th_m = 0;    % CG off set along pole
g = 9.81;

% store system parameters 
params = [Ip;
          If;
          bp;
          bf;
          m1;
          m2;
          L;
          Lm;
          g;
          th_m];

% define time span
t_final = 30;
tspan = [0 t_final];

% define intital conditions
x0 = [pi;   % theta pole
      0;     % theta flywheel
      0;     % velocity pole
      0;];  % velocity flywheel

% forward integrate
[t,x] = ode45(@(t,x)fLIP_dyna(t,x,params), tspan, x0); 

% draw fLIP
figure(1);
[row, col] = size(x);
for i = 1:row

    fLIP_draw(t(i),x(i,:),params)
end

% %% make fLIP movie
% 
% for i = 1:row
%   figure(1) ;
%   fLIP_draw(t(i),x(i,:),params)
%   F(i) = getframe(gcf) ;
%   drawnow
% end

% % create the video writer with fps
% writerObj = VideoWriter('fLIP_video.avi');
% writerObj.FrameRate = 28;
% 
% % open the video writer
% open(writerObj);
% 
% % write the frames to the video
% for i=1:length(F)
% 
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% 
% end
% 
% % close the writer object
% close(writerObj);


%% plot dynamics info
close all; clf;

% plot results
lables = ["Angle, $\theta$ [rad]";
          "Angle, $\theta_f$ [rad]";
          "Vel., $\dot{\theta}$ [rad/s]";
          "Vel., $\dot{\theta_f}$ [rad/s]"];

for i = 1:4
    subplot(2,2,i)
    plot(t, x(:,i),"Color",'blue','LineWidth',2)
    xlabel("Time, $t$ [sec]",Interpreter="latex")
    ylabel(lables(i),Interpreter="latex")
    set(gcf, 'Position', [0 0 1200 1200])
    set(gca,'TickLabelInterpreter','latex')
    grid on;
end

