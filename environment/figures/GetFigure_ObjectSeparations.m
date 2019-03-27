% GET THE AGENT SEPERATION SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = GetFigure_ObjectSeparations(SIM,DATA,currentFigure,objectNum)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

% INPUT HANDLING
if SIM.totalObjects - SIM.totalWaypoints < 2
    warning('There must be at least two collidable objects to plot seperation data.\n');
    figureHandle = [];
    return
end

% GENERATE THE SEPERATION DATA FOR A GIVEN AGENT
collidableObjectMETA = SIM.OBJECTS([SIM.OBJECTS.type] ~= OMAS_objectType.waypoint);   % The agent set
objectSubject = collidableObjectMETA(objectNum);                                      % The agent requested
objectMETA = SIM.OBJECTS([SIM.globalIDvector == objectSubject.objectID]);

agentLabel = sprintf('[ID-%s] %s',num2str(objectSubject.objectID),objectSubject.name);
figurePath = strcat(SIM.outputPath,sprintf('separations_%s ',agentLabel));

% FIGURE META DATA
figureHandle = figure('Name',agentLabel);                                   % Tab label
setappdata(figureHandle,'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'Visible','off');

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[agentStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,objectMETA.objectID,inf);

% CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
legendEntries = cell(1,(DATA.totalAgents-1));
legendCounter = 1; maxSeriesSeperation = 1;                                % Initial y-axis limit

ax = axes(figureHandle);
hold on;

% MOVE THROUGH OTHER OBJECTS
for ID2 = 1:DATA.totalObjects
    % DEFINE CONDITIONS FOR PLOTTING
    plotCondition = objectSubject.objectID ~= SIM.OBJECTS(ID2).objectID && ~strcmpi(SIM.OBJECTS(ID2).type,'WAYPOINT');
    if plotCondition
        % Generate legend entry
        legendEntries(legendCounter) = {sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(ID2).objectID),num2str(SIM.OBJECTS(ID2).name))};
        % Container for the timeseries data
        ABSseperationTimeSeries = zeros(1,size(agentStates,2));
        % Get comparative state data
        [objectStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID2).objectID,inf);
        % Reset collision instance to 0
        collisionInstance = 0;
        for simStep = 1:size(ABSseperationTimeSeries,2)
            % CALCULATE THE ABSOLUTE OBSTACLE SEPERATION TIMESERIES
            collisionCondition = objectSubject.radius + SIM.OBJECTS(ID2).radius;                     % Define the collision condition (physical seperation)
            objectSeperation = sqrt(sum((agentStates(1:3,simStep) - objectStates(1:3,simStep)).^2)); % Define object absolute (positional seperation)
            objectSeperation = objectSeperation; %- collisionCondition;                                % Subtract the physical seperation requirement from the positional seperation
            
            if 0 >= objectSeperation && collisionInstance == 0
                % OBJECT COLLISION INSTANCE
                collisionInstance = 1;
                % ADD COLLISION ANNOTATION
                annoString = strcat('Collision-',SIM.OBJECTS(ID2).name);
                xCoord = DATA.timeVector(simStep);
                yCoord = objectSeperation;
                text(ax,(xCoord),(yCoord),'|','fontsize',20,'Color',SIM.OBJECTS(ID2).colour); % '\uparrow'
                text(ax,xCoord,(yCoord + 1),annoString);
            end
            ABSseperationTimeSeries(simStep) = objectSeperation;
        end
        
        % Calculate seperation axis limit
        maxSeperation = max(ABSseperationTimeSeries); % Get the maximum seperation value
        if maxSeriesSeperation < maxSeperation
            maxSeriesSeperation = maxSeperation;
        end
        
        %% PLOT THE SEPERATION DATA ON CURRENT FIGURE
        plot(ax,DATA.timeVector,ABSseperationTimeSeries,...
            'LineStyle','-',...
            'LineWidth',DATA.figureProperties.lineWidth,...
            'Color',SIM.OBJECTS(ID2).colour);
        
        legendCounter = legendCounter + 1;                             % Increment legend entry
    end
end
% THE COLLISION CONDITION
refHandle = refline(ax,0,collisionCondition);                                % Adds a reference line with slope m and intercept b to the current axes.
set(refHandle,'color','k',...
    'LineStyle','--',...
    'LineWidth',DATA.figureProperties.lineWidth);
legendEntries{legendCounter} = 'Collision Boundary';

% ADD FINAL PLOT ATTRIBUTES
% Title 
title(ax,sprintf('Object seperations for agent %s',agentLabel),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize);
% X-axes
xlabel(ax,'t (s)',...
    'Interpreter', DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize);
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
% Y-axes
ylabel(ax,'Seperation(m)',...
    'Interpreter', DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize);
ylim(ax,[(-0.1*maxSeriesSeperation) (1.1*maxSeriesSeperation)]);
% Axes properties
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight,...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
% Legend
legend(ax,legendEntries,...
        'Interpreter',DATA.figureProperties.interpreter,...
        'Location','northeastoutside');
grid on; box on; grid minor;
% Show the timestep difference in the figure
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):SIM.TIME.dt:ax.XAxis.Limits(2);
drawnow;
hold off;

set(figureHandle,'Visible','on');

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

set(figureHandle,'Visible','off');

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;      
end