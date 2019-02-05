% GET EVENTS STATISTICS OVERVIEW PANEL
function [currentFigure,figureHandle] = get_eventOverview(SIM,DATA,currentFigure)
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
eventFields = fieldnames(DATA.eventTimeSeries); % Get the names of the event fields
for name = 1:length(eventFields)
    fieldString = char(eventFields(name));      % Convert name to string
    dataMatrix = horzcat(dataMatrix,transpose(DATA.eventTimeSeries.(fieldString)));
    % DISPLAY DATA
    eventSummary = sprintf('%s: %s',char(DATA.figureProperties.legendSet_events{name})...
                                         ,num2str(length(DATA.events.(fieldString))));
    overviewInfo = horzcat(overviewInfo,eventSummary);
end

% CREATE FIGURE HANDLE
figureHandle = figure('Name','OpenMAS Event Overview');
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.backGroundColor); % Background colour 

% MAXIMISE GRAPH SIZE IN WINDOW
% setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(gca, 'ActivePositionProperty', 'OuterPosition');
set(gca,'outerposition',[0.0 0.0 1 1])

subplot(2,3,[1 5]);
barChart = bar(DATA.timeVector(1:SIM.TIME.endStep),dataMatrix(1:SIM.TIME.endStep,:),'stacked'); % 'grouped' % Plot datamatrix against the timevector
titleStr = sprintf('Event occurance over the %ss period)',num2str(DATA.timeVector(end)));
title(titleStr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
xlabel('Time(s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
ylabel('Number of Events','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
xlim([-SIM.TIME.dt,inf])
grid on;
set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
      'fontWeight',DATA.figureProperties.fontWeight);
  
legend(DATA.figureProperties.legendSet_events); 

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
subplot(2,3,[3 6]);
axis off;
columnXOffset = 0.68;
columnWidth = 0.30;
textboxPosition = 0.6;
textBoxHeight = 0.35;
tbox = annotation(figureHandle,'textbox',[columnXOffset textboxPosition columnWidth textBoxHeight],'Units','Normalized','String',overviewInfo); 
set(tbox,'BackgroundColor','w',...
               'EdgeColor','k',...
                  'Margin', 5,...
                'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);

% CREATE THE COLLISION DATA CONTENTS
agentSummary    = sprintf('Agent Collisions: %s/%s (%3.0f%%)',num2str(DATA.collisions),num2str(DATA.totalAgents),DATA.collisionPercentage);
waypointSummary = sprintf('Waypoints: %s/%s (%3.0f%%)',num2str(DATA.waypointsAchieved),num2str(DATA.totalWaypoints),DATA.waypointPercentage);
statsInfo = {'Scenario Summary','',...
             agentSummary,waypointSummary};

% CREATE WAYPOINT DATA OVERVIEW
% tbox = annotation(figureHandle,'textbox',[0.65 0.11 .32 .4],'Units','Normalized','String',statsInfo); 
tbox = annotation(figureHandle,'textbox',[columnXOffset textboxPosition/3 columnWidth textBoxHeight],'Units','Normalized','String',statsInfo); 
set(tbox,'BackgroundColor','w',...
               'EdgeColor','k',...
                  'Margin', 5,...
                'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'eventOverview.fig'));
% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end