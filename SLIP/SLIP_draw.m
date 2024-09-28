% draw SLIP

function drawSLIP(t,x,params)
    
    % animate
    tic
    t_now = t(1);
    ind = 0;
    while t_now < t(end)
        if ind==length(t)
            break
        end
        t_now = toc;
        ind = find(t>t_now,1);
        draw(t(ind),x(ind,:),params)
    end
end

function draw(t,x,params)
    
    % position of the mass
    px = x(1);
    pz = x(2);

    % visulization setting
    wr = 0.1; % ball radius
    W = 0.1;        % base width
    H = 0.025;      % base height

    % draw pole
    plot([0,px],[0,pz],'k','LineWidth',2.5)

    % draw base
    rectangle('Position', ...
          [-W/2, -H/2, W, H], ...
          'Curvature',.4, ...
          'FaceColor', "#8c8c8c")  

    % draw ball
    rectangle('Position', ...
          [px-wr/2,pz-wr/2,wr,wr], ...
          'Curvature',1,...
          'FaceColor',"#9c1e03")

    % insert time
    time = sprintf("Time = %.2f",t);
    title(time,'Interpreter','latex')

    % set drawing boundaries
    max_radius = params.l0*1.5;
    xlim([-max_radius max_radius]); 
    ylim([-max_radius max_radius]);

    % box off
    drawnow
    set(gcf,'position',[0,0,1000,1000])
    hold off
    
end