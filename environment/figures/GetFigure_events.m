% GET EVENTS STATISTICS OVERVIEW PANEL
function [currentFigure,figureHandle] = GetFigure_events(SIM,DATA,currentFigure)
% This function generates a complete event summary plot
% INPUTS:
% SIM               - The simulation meta structure
% DATA              - The output data structure
% - eventTimeSeries - The event substructure
% currentFigure     - The current plot number

% OUTPUTS:
% currentFigure - The updated plot number 
% figureHandle  - The axes handle

% INPUT HANDLING
if ~isfield(DATA,'totalEvents')
    fprintf('[%s]\tThere was no EVENT data to plot.\n',SIM.phase);
    figureHandle = [];
    return
end

% GET EVENTS TIMESERIES DATA
dataMatrix = []; 
overviewInfo = {'Event Summary','',sprintf('Total Events: %s',num2str(DATA.totalEvents)),''};              % Plot and string data containers

% GENERATE PLOT-SPECIFIC DATA
eventFields = fieldnames(DATA.eventTimeSeries);                            % Get the names of the event fields
for name = 1:length(eventFields)
    fieldString = eventFields{name};                                       % Convert name to string
    dataMatrix = horzcat(dataMatrix,transpose(DATA.eventTimeSeries.(fieldString)));
    
    % DISPLAY DATA
    eventFields{name} = strrep(fieldString,'_',' ');                       % Remove characters for display
    eventSummary = sprintf('%s: %s',eventFields{name},num2str(length(DATA.events.(fieldString))));
    overviewInfo = horzcat(overviewInfo,eventSummary);
end

% FIGURE META PROPERTIES
figurePath = strcat(SIM.outputPath,'event-overview');
figureHandle = figure('Name','OpenMAS event overview');
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

% GENERATE THE SUB PLOTS (SPECIFIC ORIENTATIONS)
subHandles(1) = subplot(2,4,[1 7]);
subHandles(2) = subplot(2,4,4);
subHandles(3) = subplot(2,4,8); 

% EVENT OCCURANCE BAR CHART
h = subplot(subHandles(1));
barChart = bar(h,DATA.timeVector(1:SIM.TIME.endStep),dataMatrix(1:SIM.TIME.endStep,:),'stacked'); % 'grouped' % Plot datamatrix against the timevector

% Title
title(h,...
    sprintf('Event occurance over a period of %ss',num2str(SIM.TIME.endTime)),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-Label
xlabel(h,...
    't (s)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
xlim([-SIM.TIME.dt,SIM.TIME.endTime])
% Y-Label
ylabel(h,...
    'Number of Events',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Axes
set(h,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor);
% Legend
legend(h,eventFields,...
    'Location','northeastoutside',...
    'fontname',DATA.figureProperties.fontName,...
    'Interpreter',DATA.figureProperties.interpreter);
grid on; grid minor; box on;

% ASSIGN SPECIFIC COLOURS TO KEY PROPERTIES
colourVector = [0 0 0;     % 'event' type
                0 0 1;     % 'detection' type
                1 0.5 0;   % 'warning' type
                1 0 0;     % 'collision' type
                0 1 0;     % 'waypoint' type
                0 1 1;     % 'null_detection' type
                1 1 0;     % 'null_warning' type  
                1 0 1;     % 'null_collision' type
                0.2 1 0.2];% 'null_waypoint' type   
            
% ALLOCATE COLOURS TO THE FIGURE BARS
[eventEnum,~] = enumeration('eventType');
for i = 1:numel(eventEnum)
    enumIndex = double(eventEnum(i));
    barChart(enumIndex + 1).FaceColor = colourVector(i,:);
end            
            
% BUILD THE NEIGHBOURING TEXT BLOCKS
h = subplot(subHandles(2));   
axis off;
columnXOffset = 0.72;
columnWidth = 0.250;
textboxYOffset = 0.42;
textBoxHeight = 0.5;
tboxA = annotation(figureHandle,'textbox',...
    [columnXOffset textboxYOffset columnWidth textBoxHeight],...
    'Units','Normalized',...
    'String',overviewInfo);
set(tboxA,...
    'Interpreter',DATA.figureProperties.interpreter,...
    'BackgroundColor','w',...
    'EdgeColor','k',...
    'Margin',5,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight);

% BUILD THE NEIGHBOURING TEXT BLOCKS
h = subplot(subHandles(3));
axis off;  
% CREATE THE COLLISION DATA CONTENTS
agentSummary    = sprintf('Agent Collisions: %d/%d',DATA.collisions,DATA.totalAgents); % DATA.collisionPercentage
waypointSummary = sprintf('Waypoints: %d/%d',DATA.waypointsAchieved,DATA.totalWaypoints); % DATA.waypointPercentage
statsString     = {'Scenario Summary','',agentSummary,waypointSummary};

% CREATE WAYPOINT DATA OVERVIEW
textboxYOffset = 0.11;
textBoxHeight = 0.3;
tboxB = annotation(figureHandle,'textbox',...
    [columnXOffset textboxYOffset columnWidth textBoxHeight],...
    'Units','Normalized',...
    'String',statsString);
set(tboxB,...
    'Interpreter',DATA.figureProperties.interpreter,...
    'BackgroundColor','w',...
    'EdgeColor','k',...
    'Margin',5,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);  

% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end
