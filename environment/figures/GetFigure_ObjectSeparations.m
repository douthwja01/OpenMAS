%% GET THE AGENT SEPERATION SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = GetFigure_objectSeparations(SIM,DATA,currentFigure,objectNum)
% Draws the separations between the 'subject object', its associated waypoints 
% and obstacles.

% Input sanity check
if sum([SIM.OBJECTS.type] ~= OMAS_objectType.waypoint) < 2
    warning('There must be at least two collidable objects to plot separation data.');
    figureHandle = [];
    return
end
% Subject object reference
subjectMETA = SIM.OBJECTS(objectNum);   % META data associated with the 'subject' object
% Second object references
collidableMETA = SIM.OBJECTS([SIM.OBJECTS.hitBox] ~= OMAS_hitBoxType.none);         % Only collidable objects
collidableMETA = collidableMETA([collidableMETA.objectID] ~= subjectMETA.objectID); % Not the same object 
collidableMETA = collidableMETA([collidableMETA.type] ~= OMAS_objectType.waypoint); % That are not waypoints

% Agent label
agentLabel = sprintf('[ID-%d] %s',subjectMETA.objectID,subjectMETA.name);

% Figure META data
figurePath = strcat(SIM.outputPath,sprintf('separations_%s ',agentLabel));
figureHandle = figure('Name',agentLabel);                                   % Tab label
setappdata(figureHandle,'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'Visible','off');
ax = axes(figureHandle);
hold on;

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[subjectStates] = OMAS_getTrajectoryData_mex(...
    DATA.globalTrajectories,...
    SIM.globalIDvector,...
    subjectMETA.objectID,...
    inf);

% Constants
runSteps = size(subjectStates,2);
legendEntries = cell(numel(collidableMETA),1);
maxSeriesSeperation = 1;  

% All objects in the collidable set are to be plotted
for i = 1:numel(collidableMETA)
    % Modify the legend entries
    legendEntries{i} = sprintf('[ID-%d] %s',collidableMETA(i).objectID,collidableMETA(i).name);
    % Container for the timeseries data
    ABSseperationTimeSeries = zeros(1,runSteps);
    % Get comparative state data
    collidableStates = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,collidableMETA(i).objectID,inf);
    % Reset collision instance to 0
    collisionInstance = 0;
    for simStep = 1:runSteps
        % Define the collision condition
        collisionCondition = subjectMETA.radius + collidableMETA(i).radius;                     % Define the collision condition (physical seperation)
        % Calculate the seperatation at each timestep
        objectSeperation = sqrt(sum((subjectStates(1:3,simStep) - collidableStates(1:3,simStep)).^2)); % Define object absolute (positional seperation)
        
        if 0 >= objectSeperation && collisionInstance == 0
            % Object collision instance
            collisionInstance = 1;
            % Add collision annotation
            xCoord = DATA.timeVector(simStep);
            yCoord = objectSeperation;
            text(ax,xCoord,yCoord,'|',...
                'fontname',DATA.figureProperties.fontName,...,...
                'fontsize',20,...
                'Color',collidableMETA(i).colour); % '\uparrow'
            text(ax,xCoord,(yCoord + 1),strcat('Collision-',collidableMETA(i).name));
        end
        ABSseperationTimeSeries(simStep) = objectSeperation;
    end
    
    % Calculate seperation axis limit
    maxSeparation = max(ABSseperationTimeSeries); % Get the maximum seperation value
    if maxSeriesSeperation < maxSeparation
        maxSeriesSeperation = maxSeparation;
    end
    
    % Plot the separation data
    l = plot(ax,DATA.timeVector(1:SIM.TIME.endStep),ABSseperationTimeSeries(1:SIM.TIME.endStep));
    set(l,'LineWidth',DATA.figureProperties.lineWidth);
    set(l,'Color',collidableMETA(i).colour);
    % Data presentation
    switch collidableMETA(i).type
        case OMAS_objectType.waypoint
            set(l,'LineStyle','--');
        otherwise
            set(l,'LineStyle','-');
    end
end

%% Additional plot refinements
% The collision condition
refHandle = refline(ax,0,collisionCondition);                                % Adds a reference line with slope m and intercept b to the current axes.
set(refHandle,'Color','k');
set(refHandle,'LineStyle','--');
set(refHandle,'LineWidth',DATA.figureProperties.lineWidth/2);
legendEntries = vertcat(legendEntries,'Collision Boundary');
% Title 
title(ax,sprintf('Object separations for agent %s',agentLabel),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize);
% X-axes
xlabel(ax,'t (s)',...
    'Interpreter', DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize);
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
% Y-axes
ylabel(ax,'Separation (m)',...
    'Interpreter', DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize);
ylim(ax,[0 (1.1*maxSeriesSeperation)]);
% Axes properties
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight,...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
% Legend
legend(ax,legendEntries,...
       'Interpreter',DATA.figureProperties.interpreter,...
       'fontname',DATA.figureProperties.fontName,...
       'Location','northeast');
grid on; box on; grid minor;
% Show the timestep difference in the figure
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):SIM.TIME.dt:ax.XAxis.Limits(2);
drawnow;
hold off;
set(figureHandle,'Visible','on');

% Save the .fig by default
savefig(figureHandle,figurePath);      

% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end

set(figureHandle,'Visible','off');

% Modify offset if required
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
% Incremente the figure index
currentFigure = currentFigure + 1;      
end