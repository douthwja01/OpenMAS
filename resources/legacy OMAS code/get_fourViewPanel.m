% GET THE 3 PLAN + ISOMETRIC ANIMATED TRAILS FIGURE
function [currentFigure,figureHandle] = get_fourViewPanel(SIM,DATA,currentFigure)
% This function generates a four-view animated panel showing the motion of
% all objects in the simulation.

% BUILD FIGURE
% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Four Viewpoint Panel');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]

% GET THE LAST INSTANCE OF A VALID STATE
for step = 1:SIM.TIME.endStep
    try
        % CLEAR MAIN PLOT
        cla reset;
        % CLEAR SUBPLOTS
        cla(s1);cla(s2);cla(s3);cla(s4);
    catch
    end
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object
        activeStates = size(objectStates,2);
        if step > activeStates
            tailLength = DATA.figureProperties.tailLength/SIM.TIME.dt;
            if activeStates < tailLength
                tailLength = activeStates;
            end
           % OBJECT HAS BEEN REMOVED
           markerPosition = objectStates(1:3,end);
           tailTrajectories = objectStates(1:3,((end+1)-tailLength):end);
        elseif step > (DATA.figureProperties.tailLength/SIM.TIME.dt)
           % ONLY THE REQUESTED TAIL LENGTH
           markerPosition = objectStates(1:3,step);
           tailTrajectories = objectStates(1:3,(step - DATA.figureProperties.tailLength/SIM.TIME.dt):step);
        else
            % TRANSITIONING PERIOD
            markerPosition = objectStates(1:3,step);
            tailTrajectories = objectStates(1:3,1:step);
        end
        
        setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
        % PLOT ONE (ISOMETRIC) ////////////////////////////////////////////
        s1 = subplot(4,4,[1 6]);
        hold on; grid on; box on;
        % PLOT THE HEAD POSITIONS AT STEP VALUE
        plot3(s1,markerPosition(1),markerPosition(2),markerPosition(3),...
            'Marker',SIM.OBJECTS(indexValue).symbol,...
            'MarkerSize',DATA.figureProperties.MarkerSize,...
            'MarkerFaceColor',SIM.OBJECTS(indexValue).colour,...
            'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
            'Color',SIM.OBJECTS(indexValue).colour);
        % PLOT THE TAIL DATA
        plot3(s1,tailTrajectories(1,:),tailTrajectories(2,:),tailTrajectories(3,:),...
            'LineStyle',DATA.figureProperties.LineStyle,...
            'LineWidth',DATA.figureProperties.LineWidth,...
            'Color',SIM.OBJECTS(indexValue).colour);
        title('Isometric View','fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);        
        ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
        zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        view([-45 50]);
        hold off;
        
        % PLOT TWO ( PLANE XY ) ///////////////////////////////////////////
        s2 = subplot(4,4,[3 8]);  % Plot each state graph
        hold on; grid on; box on;
        % PLOT THE HEAD POSITIONS AT STEP VALUE
        plot(s2,markerPosition(1),markerPosition(2),...
            'Marker',SIM.OBJECTS(indexValue).symbol,...
            'MarkerSize',DATA.figureProperties.MarkerSize,...
            'MarkerFaceColor',SIM.OBJECTS(indexValue).colour,...
            'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
            'Color',SIM.OBJECTS(indexValue).colour);
        % PLOT THE TAIL DATA
        plot(s2,tailTrajectories(1,:),tailTrajectories(2,:),...
            'LineStyle',DATA.figureProperties.LineStyle,...
            'LineWidth',DATA.figureProperties.LineWidth,...
            'Color',SIM.OBJECTS(indexValue).colour);
        title('XY Plane','fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        hold off;
        legend('off')
        
        % PLOT THREE ( PLANE YZ ) /////////////////////////////////////////
        s3 = subplot(4,4,[9 14]);  % Plot each state graph
        hold on; grid on; box on;
        % PLOT THE HEAD POSITIONS AT STEP VALUE
        plot(s3,markerPosition(2),markerPosition(3),...
            'Marker',SIM.OBJECTS(indexValue).symbol,...
            'MarkerSize',DATA.figureProperties.MarkerSize,...
            'MarkerFaceColor',SIM.OBJECTS(indexValue).colour,...
            'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
            'Color',SIM.OBJECTS(indexValue).colour);
        % PLOT THE TAIL DATA
        plot(s3,tailTrajectories(2,:),tailTrajectories(3,:),...
            'LineStyle',DATA.figureProperties.LineStyle,...
            'LineWidth',DATA.figureProperties.LineWidth,...
            'Color',SIM.OBJECTS(indexValue).colour);
        title('YZ Plane','fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        hold off;
        legend('off')
        
        % PLOT FOUR ( PLANE XZ ) //////////////////////////////////////////
        s4 = subplot(4,4,[11 16]);  % Plot each state graph
        hold on; grid on; box on;
        % PLOT THE HEAD POSITIONS AT STEP VALUE
        plot(s4,markerPosition(1),markerPosition(3),...
            'Marker',SIM.OBJECTS(indexValue).symbol,...
            'MarkerSize',DATA.figureProperties.MarkerSize,...
            'MarkerFaceColor',SIM.OBJECTS(indexValue).colour,...
            'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
            'Color',SIM.OBJECTS(indexValue).colour);
        % PLOT THE TAIL DATA
        plot(s4,tailTrajectories(1,:),tailTrajectories(3,:),...
            'LineStyle',DATA.figureProperties.LineStyle,...
            'LineWidth',DATA.figureProperties.LineWidth,...
            'Color',SIM.OBJECTS(indexValue).colour);
        title('XZ Plane','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.titleFontSize);
        xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
        set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        hold off;
        legend('off')
    end
    
    % DRAW AND COLLATE GIF FRAME
    frame = getframe(figureHandle);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % APPEND THE FRAMES
    if step == 1
        imwrite(imind,cm,strcat(SIM.outputPath,'fourView','.gif')...
            ,'gif', 'Loopcount',inf,'DelayTime',SIM.TIME.dt);
    else
        imwrite(imind,cm,strcat(SIM.outputPath,'fourView','.gif')...
            ,'gif','WriteMode','append','DelayTime',SIM.TIME.dt); 
    end
end
% legend(s1,DATA.figureProperties.legendEntries,'location','BestOutside');

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end