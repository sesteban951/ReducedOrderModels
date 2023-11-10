%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fLIP DRAW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw(t,x,params)

% system parameters
L = params(7);
Lm = params(8);
th_m = params(9);

% get state data
th = x(:,1);
th_f = x(:,2);          % position of the flywheel doesnt affect dynamics

% % fix for plotting
% th = th;   
% th_m = th_m; 

% flip beacuse sign convention is wrong for some reason

% flywheel and base dimensions
wr = 0.3 * L;     % flywheel radius
cgr = 0.4 * wr;   % CG marker radius
W = 0.4  * L;     % base width
H = 0.025 * L;    % base height

% position of flywheel
px = L * sin(th);    
py = L * cos(th);

% position of CG
px_cg = Lm * sin(th + th_m);
py_cg = Lm * cos(th + th_m);

% draw fLIP system
% hold on;

% draw pole
plot([0,px],[0,py],'w','LineWidth',2.5)

% draw base
rectangle('Position', ...
          [-W/2, -H/2, W, H], ...
          'Curvature',.4, ...
          'FaceColor', "#8c8c8c")   

% draw flywheel
rectangle('Position', ...
          [px-wr/2,py-wr/2,wr,wr], ...
          'Curvature',1,...
          'FaceColor',"#9c1e03")   

% draw CG
rectangle('Position', ...
          [px_cg-cgr/2,py_cg-cgr/2,cgr,cgr], ...
          'Curvature',1,...
          'FaceColor',"#47ed0a")  

% draw flywheel position
line([px,px+wr*cos(th_f)/2],[py,py+wr*sin(th_f)/2], ...
    'Color','w', ...
    'LineWidth',2)

% insert time
text(0,-0.1,"t = " + string(t),'Color','w','FontSize',20)

% set drawing boundaries
max_radius = (L + wr)*1.1;
xlim([-max_radius max_radius]); 
ylim([-max_radius max_radius]);
set(gca,'Color','k','XColor','w','YColor','w')
set(gcf,'Color','k')
set(gcf,'InvertHardcopy','off')   
axis equal

% box off
drawnow
set(gcf, 'Position', [0 0 1000 1000])
hold off
