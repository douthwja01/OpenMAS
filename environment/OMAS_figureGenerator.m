%% OPENMAS FIGURE GENERATOR (OMAS_figureGenerator.m) %%%%%%%%%%%%%%%%
% This script contains an index of the figures that can be requested from
% the simulation using the SIM.figures attribute.

% Author: James A. Douthwaite 10/10/2016

%% MASTER FUNCTION
function [figureNumber] = OMAS_figureGenerator(SIM,DATA,figureNumber,figureLabel)
% INPUTS:
% DATA         - The simulation output DATA structure
% figureNumber - The current figure number
% figureFlag   - The figure identifier
% OUTPUT:
% figureNumber - The updated figure number
   
% DETERMINE WHICH FIGURE IS TO BE GENERATED
switch upper(char(figureLabel))
    case 'ALL'
        fprintf('[%s]\tAll output figures requested.\n',SIM.phase);
        % MOVE THROUGH THE COMPLETE FIGURE VECTOR
        figureVector = {'EVENTS','COLLISIONS','TRAJECTORIES','SEPERATIONS',...
                        'INPUTS','ISOMETRIC','4VIEW','GIF','VIDEO','TIMES'};
        for fig = 1:length(figureVector)
            [figureNumber] = OMAS_figureGenerator(SIM,DATA,figureNumber,figureVector{fig});
        end
        close all;
    case 'EVENTS'
        fprintf('[%s]\tGenerating the event overview figure.\n',SIM.phase);
        [figureNumber,~] = get_eventOverview(SIM,DATA,figureNumber);    
    case 'COLLISIONS'
        fprintf('[%s]\tGenerating collision figure set.\n',SIM.phase);
        [figureNumber,~] = get_collisionFigures(SIM,DATA,figureNumber);
    case 'TRAJECTORIES'
        fprintf('[%s]\tGenerating global trajectory figure.\n',SIM.phase);
        figureSet = [];
        for objectID = 1:DATA.totalObjects
            [figureNumber,figureSet(objectID)] = get_objectTrajectory(SIM,DATA,figureNumber,objectID);
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = OMAS_figureTabUtility(figureSet,'OpenMAS Trajectory Overview');
        % MAXIMISE GRAPH SIZE IN WINDOW
%         setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
        set(windowHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-80)]);
        % SAVE THE OUTPUT FIGURE
        savefig(windowHandle,strcat(SIM.outputPath,'trajectoryFigure.fig'));
    case 'SEPERATIONS' 
        fprintf('[%s]\tGenerating object trajectory seperation figure.\n',SIM.phase);
        [figureNumber,~] = get_agentSeperationFigure(SIM,DATA,figureNumber);
    case 'INPUTS'
        fprintf('[%s}\tGenerating control input figure.\n',SIM.phase);
        [figureNumber,~] = get_agentControlInputs(SIM,DATA,figureNumber);
    case 'ISOMETRIC' 
        fprintf('[%s]\tGenerating isometric trajectory figure.\n',SIM.phase);
        [figureNumber,~] = get_isometricFigure(SIM,DATA,figureNumber);
    case '4VIEW'
        fprintf('[%s]\tGenerating four-view panel & gif file.\n',SIM.phase);
        [figureNumber,~] = get_fourViewPanel(SIM,DATA,figureNumber);
    case 'GIF'
        fprintf('[%s]\tGenerating trajectory gif file.\n',SIM.phase);
        [figureNumber,~] = get_isometricGif(SIM,DATA,figureNumber);
    case 'VIDEO'
        fprintf('[%s]\tGenerating trajectory avi file.\n',SIM.phase);
        [figureNumber,~] = get_isometricVideo(SIM,DATA,figureNumber);
    case 'TIMES'
        fprintf('[%s]\tGenerating computation timeseries figure.\n',SIM.phase);
        [figureNumber,~] = get_computationTimes(SIM,DATA,figureNumber);    
    otherwise
        warningStr = sprintf('[ERROR] Did not recognise figure request: "%s"',char(figureLabel));
        warning(warningStr);
end
end

%% FIGURE GENERATION FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COLLISION DATA FIGURE GENERATION
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
    eventSummary = sprintf('%s: %s',char(DATA.figureProperties.legendSet_events{name})...
                                         ,num2str(length(DATA.events.(fieldString))));
    overviewInfo = horzcat(overviewInfo,eventSummary);
end

% CREATE FIGURE HANDLE
figureHandle = figure('Name','OpenMAS Event Overview');
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

subplot(2,3,[1 5]);
barChart = bar(DATA.timeVector,dataMatrix,'stacked'); % 'grouped' % Plot datamatrix against the timevector
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
tbox = annotation(figureHandle,'textbox',[columnXOffset 0.62 columnWidth .30],'Units','Normalized','String',overviewInfo); 
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
tbox = annotation(figureHandle,'textbox',[columnXOffset 0.4 columnWidth .20],'Units','Normalized','String',statsInfo); 
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
% PLOT THE COLLISION EVENT DESCRIPTION FIGURES
function [currentFigure,tabbedFigureHandle] = get_collisionFigures(SIM,DATA,currentFigure)
% This function is designed to move through the collsion event history and
% generate figure for each of the collisions.

% CHECK COLLISIONS OCCURED
if ~isfield(DATA,'events') || ~isfield(DATA.events,'collisions') 
    fprintf('[%s]\t...No collision data available.\n',SIM.phase);
    tabbedFigureHandle = 0;
    return
end

figureSet = [];
for collisionNumber = 1:size(DATA.events.collisions,1)
    collisionEvent = DATA.events.collisions(collisionNumber);              % Get the collision event from the event history structure.
    collisionTime = collisionEvent.time;                                   % Define the time of the collision
    
    % Generate the figure
    tabStr = sprintf('Collision %s @ %s',num2str(collisionNumber),num2str(collisionTime));
    newFigureHandle = figure('Name',tabStr);
    set(newFigureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
    
    hold on; grid on; box on;                                              % Assign appropriate figure properties                                         
    axis equal;
    view([70 25]); 
    titlestr = sprintf('Collision Event between %s and %s at t=%s',...     % Declare title string for figure
                       collisionEvent.name_A,collisionEvent.name_B,num2str(collisionEvent.time));
    % Get trajectory data
    [objectAStates] = OMAS_getTrajectoryData(DATA,collisionEvent.objectID_A);              % Get the first object trajectory data
    [objectBStates] = OMAS_getTrajectoryData(DATA,collisionEvent.objectID_B);              % Get the second object trajectory data
    
    %% PLOT THE TRAILS
    plot3(objectAStates(1,:),objectAStates(2,:),objectAStates(3,:),...
                'LineStyle',DATA.figureProperties.LineStyle,...
                'LineWidth',DATA.figureProperties.LineWidth,...
                'Color',SIM.OBJECTS(collisionEvent.objectID_A).colour);
    plot3(objectBStates(1,:),objectBStates(2,:),objectBStates(3,:),...
                'LineStyle',DATA.figureProperties.LineStyle,...
                'LineWidth',DATA.figureProperties.LineWidth,...
                'Color',SIM.OBJECTS(collisionEvent.objectID_B).colour);
     
    %% DEFINE THE COLLISION SPHERES 
    % Build object A volume representation 
    radius = SIM.OBJECTS(collisionEvent.objectID_A).size;
    [X,Y,Z] = sphere(40);
    Xa = X.*radius + collisionEvent.state_A(1);
    Ya = Y.*radius + collisionEvent.state_A(2);
    Za = Z.*radius + collisionEvent.state_A(3);    
    mesh(Xa,Ya,Za,...
        'FaceColor',SIM.OBJECTS(collisionEvent.objectID_A).colour,...
        'FaceAlpha',0.1,...
        'LineWidth',DATA.figureProperties.LineWidth,...
        'EdgeAlpha',0.2,...
        'edgecolor',DATA.figureProperties.MarkerEdgeColor);  % Obstacle

    % Build object B volume representation
    radius = SIM.OBJECTS(collisionEvent.objectID_B).size;
    Xb = X.*radius + collisionEvent.state_B(1);
    Yb = Y.*radius + collisionEvent.state_B(2);
    Zb = Z.*radius + collisionEvent.state_B(3);
    mesh(Xb,Yb,Zb,...
        'FaceColor',SIM.OBJECTS(collisionEvent.objectID_B).colour,...
        'FaceAlpha',0.1,...
        'LineWidth',DATA.figureProperties.LineWidth,...
        'EdgeAlpha',0.2,...
        'edgecolor',DATA.figureProperties.MarkerEdgeColor);  % Obstacle

    %% PLOT THE OBJECT VELOCITIES
    q = quiver3(collisionEvent.state_A(1),...
                collisionEvent.state_A(2),...
                collisionEvent.state_A(3),...
                collisionEvent.state_A(7),...
                collisionEvent.state_A(8),...
                collisionEvent.state_A(9),'r');
    q.AutoScaleFactor = 1;
    q.LineWidth = DATA.figureProperties.LineWidth;
    q = quiver3(collisionEvent.state_B(1),...
                collisionEvent.state_B(2),...
                collisionEvent.state_B(3),...
                collisionEvent.state_B(7),...
                collisionEvent.state_B(8),...
                collisionEvent.state_B(9),'g');
    q.AutoScaleFactor = 1;
    q.LineWidth = DATA.figureProperties.LineWidth;   
    
    % Other plot attributes
    title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
    xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
          'fontWeight',DATA.figureProperties.fontWeight);
    hold off;
    
    figureSet = vertcat(figureSet,newFigureHandle);
    % SAVE THE OUTPUT FIGURE
%     fileName = sprintf('Collision(%s)[%s-%s].fig',num2str(collisionNumber),...
%                         collisionEvent.name_A,collisionEvent.name_B);
%     savefig(figureHandle,strcat(SIM.outputPath,fileName));
    
    % Iterate the figure counter for next collision instance
%     DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
%     currentFigure = currentFigure + 1;
end

% ASSEMBLE TABBED FIGURE
tabbedFigureHandle = OMAS_figureTabUtility(figureSet,'OpenMAS Collision Overview');

% SAVE THE OUTPUT FIGURE
savefig(tabbedFigureHandle,strcat(SIM.outputPath,'collisionOverview.fig'));

currentFigure = currentFigure + 1;
end

%% TRAJECTORY DATA FIGURE GENERATION
% GET THE GLOBAL TRAJECTORY DATA FOR AN INDIVIDUAL OBJECT
function [currentFigure,figureHandle] = get_objectTrajectory(SIM,DATA,currentFigure,objectNum)
% Declare title string for figure  
titlestr = sprintf('Global trajectory data for %s[ID:%s] over a period of %ss'...
                    ,SIM.OBJECTS(objectNum).name,num2str(objectNum),num2str(DATA.timeVector(end)));

% CONFIGURE THE PLOT ATTRIBUTES
% figureHandle = figure(currentFigure);
figureHandle = figure('Name',SIM.OBJECTS(objectNum).name);
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]);
plotCellWidth = 4; plotCellA = 1;                                          % The width of each figure, the start of the plot

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[objectStates] = OMAS_getTrajectoryData(DATA,objectNum);              % Get the object data
% STATE NAME VECTOR
stateTags = {'x_{m}','y_{m}','z_{m}',...
             'u_{m/s}','v_{m/s}','w_{m/s}',...
             'q0','q1','q2','q3'};

setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08,0.08,0.90,0.90]);         

% FOR EACH OBJECTS PERSPECTIVE
stateVectorLength = size(objectStates,1);
for stateNum = 1:stateVectorLength
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    stateTrajectory = objectStates(stateNum,:);
    % GET THE Y-AXIS STATE LABEL
    yString = stateTags{stateNum};
        
    % CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
    plotCellB = stateNum*plotCellWidth;                                    % The end of the plot
    plotLocation = subplot(stateVectorLength,plotCellWidth,[plotCellA plotCellB]);
    hold on;
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);                                                   % Append title to first subplot
    end
    
    %% PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,DATA.timeVector,stateTrajectory,...
         'LineStyle','-',...
         'LineWidth',DATA.figureProperties.LineWidth,...
         'Color','b');
    ylabel(yString,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    
    % Prevent overlap of x-label
    if plotCellB == (stateVectorLength*plotCellWidth) 
        xlabel('t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    end
    set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
          'fontWeight',DATA.figureProperties.fontWeight);
    grid on; box on;
    
    % Move to next subplot location 
    plotCellA = plotCellA + plotCellWidth;
end
hold off;
currentFigure = currentFigure + 1;
end
% GET THE OBJECT SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = get_agentSeperationFigure(SIM,DATA,currentFigure)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

if DATA.totalAgents <= 1
	warning('Single agent simulation - no seperation data available.');
    figureHandle = [];
    return
end

% Declare title string for figure  
titlestr = sprintf('Agent seperations over a period of %ss',num2str(DATA.timeVector(end)));

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Seperation Overview');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
plotCellWidth = 4; plotCellA = 1;                                          % The width of each figure, the start of the plot

% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
% FOR EACH AGENT'S PERSPECTIVE
% for ID1 = 1:DATA.totalObjects
for ID1 = 1:DATA.totalAgents
    % CHECK IF EVALUATION OBJECT IS A WAYPOINT; OMIT IF NECESSARY
    if strcmpi(SIM.OBJECTS(ID1).type,'WAYPOINT')
        continue
    end
    
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    [objectID1States] = OMAS_getTrajectoryData(DATA,ID1);              % Trajectory data per objectID
    objectID1 = SIM.OBJECTS(ID1).objectID;
    objectID1Name = SIM.OBJECTS(ID1).name;
    ystring = sprintf('%s [ID:%s] \n Seperations (m)',objectID1Name,num2str(objectID1));
    
    % CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
    legendEntries = cell(1,(DATA.totalAgents-1));
    legendCounter = 1; maxSeriesSeperation = 1;                            % Initial y-axis limit
    plotCellB = ID1*plotCellWidth;                                         % The end of the plot
    plotLocation = subplot(DATA.totalAgents,plotCellWidth,[plotCellA plotCellB]);
    hold on;
    
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);                                           % Append title to first subplot
    end
    
    % MOVE THROUGH OTHER OBJECTS
    for ID2 = 1:DATA.totalObjects
        % DEFINE CONDITIONS FOR PLOTTING
        plotCondition = SIM.OBJECTS(ID1).objectID ~= SIM.OBJECTS(ID2).objectID && ~strcmpi(SIM.OBJECTS(ID2).type,'WAYPOINT');
        if plotCondition
            objectID2 = SIM.OBJECTS(ID2).objectID;                    % Assemble the object legend information
            objectID2Name = num2str(SIM.OBJECTS(ID2).name);
            % Generate legend entry
            legendEntries(legendCounter) = {sprintf('[ID:%s] %s',num2str(objectID2),objectID2Name)};
            % Container for the timeseries data
            ABSseperationTimeSeries = zeros(1,size(objectID1States,2));
            % Get comparative state data
            [objectID2States] = OMAS_getTrajectoryData(DATA,ID2);
            
            % Reset collision instance to 0
            collisionInstance = 0;
            for simStep = 1:size(ABSseperationTimeSeries,2)
                % CALCULATE THE ABSOLUTE OBSTACLE SEPERATION TIMESERIES
                collisionCondition = SIM.OBJECTS(ID1).size + SIM.OBJECTS(ID2).size;                             % Define the collision condition (physical seperation)                      
                objectSeperation = sqrt(sum((objectID1States(1:3,simStep) - objectID2States(1:3,simStep)).^2)); % Define object absolute (positional seperation)  
                objectSeperation = objectSeperation - collisionCondition;                                       % Subtract the physical seperation requirement from the positional seperation
                
                if 0 >= objectSeperation && collisionInstance == 0
                    % OBJECT COLLISION INSTANCE
                    collisionInstance = 1;
                    % ADD COLLISION ANNOTATION
                    annoString = strcat('Collision-',objectID2Name);
                    xCoord = DATA.timeVector(simStep);
                    yCoord = objectSeperation;
                    text(plotLocation,(xCoord),(yCoord),'|','fontsize',20,'Color',SIM.OBJECTS(ID2).colour); % '\uparrow'
                    text(plotLocation,xCoord,(yCoord + 1),annoString);
                end
                ABSseperationTimeSeries(simStep) = objectSeperation;
            end
            
            % Calculate seperation axis limit
            maxSeperation = max(ABSseperationTimeSeries); % Get the maximum seperation value
            if maxSeriesSeperation < maxSeperation 
                maxSeriesSeperation = maxSeperation;
            end
            
            %% PLOT THE SEPERATION DATA ON CURRENT FIGURE
            plot(plotLocation,DATA.timeVector,ABSseperationTimeSeries,...
                        'LineStyle','-',...
                        'LineWidth',DATA.figureProperties.LineWidth,...                
                        'Color',SIM.OBJECTS(ID2).colour);
            xlabel('t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
            ylabel(ystring,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
            ylim([(-0.1*maxSeriesSeperation) (1.1*maxSeriesSeperation)]);
            grid on; box on;
            set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
                  'fontWeight',DATA.figureProperties.fontWeight);
            legendCounter = legendCounter + 1;                             % Increment legend entry
        end
    end
    
    % ADD LEGEND
    legend(legendEntries,'Location','southeast');
           
    % ADD ANNOTATION THE FOR THE VISUAL HORIZON
%     line([SIM.TIME.startTime SIM.TIME.endTime],...
%          [SIM.OBJECTS(ID1).detectionRange SIM.OBJECTS(ID1).detectionRange],...
%          'LineStyle','--','color','bl');
%     text((SIM.TIME.startTime + SIM.TIME.dt),(SIM.OBJECTS(ID1).detectionRange - 5),...
%          'Visual horizon','color','bl');   % Annotate line as visual horizon

    plotCellA = plotCellA + plotCellWidth; % Move to next subplot location
end
hold off;
% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'seperations.fig'));            % As matlab figure 
saveas(figureHandle,strcat(SIM.outputPath,'seperations'),'epsc');          % As postscript
% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;      
end
% GET THE GLOBAL SEPERATIONS OF A SINGLE AGENT/OBJECT
function [currentFigure,figureHandle] = get_objectSeperationFigure(SIM,DATA,currentFigure,objectNum)
% Declare title string for figure  
titlestr = sprintf('Global seperations of object %s over a period of %ss',...
                    SIM.OBJECTS(objectNum).name,num2str(DATA.timeVector(end)));

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name',SIM.OBJECTS(objectNum).name);
axesHandle = axes(figureHandle);
legend(axesHandle,'off')
title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize); 

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[objectAStates] = OMAS_getTrajectoryData(DATA,objectNum);            % Get the object data
objectIDA = SIM.OBJECTS(objectNum).objectID;
objectNameA = SIM.OBJECTS(objectNum).name;
yAxisLabel = {strcat(objectNameA,' [ID:',num2str(objectIDA),']'),'Seperations (m)'};   % Define y axis label

% CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
legendEntries = cell(1,(DATA.totalObjects-1));
legendCounter = 1; maxSeriesSeperation = 1;                                % Initial y-axis limit
hold on;

% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(gcf, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

for ID2 = 1:DATA.totalObjects
    if objectIDA ~= SIM.OBJECTS(ID2).objectID
        % GENERAL INFORMATION
        objectIDB = SIM.OBJECTS(ID2).objectID;                             % Assemble the object legend information
        objectNameB = num2str(SIM.OBJECTS(ID2).name);
        legendEntries(legendCounter) = {strcat('[ID:',num2str(objectIDB),'] ',objectNameB)};
        legendCounter = legendCounter + 1;                                 % Increment legend entry
        % Container for the timeseries data
        ABSseperationTimeSeries = zeros(1,size(objectAStates,2));
        % Get comparative state data
        [objectBStates] = OMAS_getTrajectoryData(DATA,ID2);
        % Reset collision instance to 0
        collisionInstance = 0;
        for simStep = 1:size(ABSseperationTimeSeries,2)
            % CALCULATE THE ABSOLUTE OBSTACLE SEPERATION TIMESERIES
            collisionCondition = SIM.OBJECTS(objectNum).size + SIM.OBJECTS(ID2).size;                   % Define the collision condition (physical seperation)
            objectSeperation = sqrt(sum((objectAStates(1:3,simStep) - objectBStates(1:3,simStep)).^2)); % Define object absolute (positional seperation)
            objectSeperation = objectSeperation - collisionCondition;                                   % Subtract the physical seperation requirement from the positional seperation
            if 0 >= objectSeperation && collisionInstance == 0
                % OBJECT COLLISION INSTANCE
                collisionInstance = 1;
                % ADD COLLISION ANNOTATION
                annoString = strcat('Collision-',objectNameB);
                xCoord = DATA.timeVector(simStep);
                yCoord = objectSeperation;
                text((xCoord),(yCoord),'|','fontsize',20,'Color',SIM.OBJECTS(ID2).colour); % '\uparrow'
                text(xCoord,(yCoord + 1),annoString);
            end
            ABSseperationTimeSeries(simStep) = objectSeperation;
        end

        % Calculate seperation axis limit
        maxSeperation = max(ABSseperationTimeSeries); % Get the maximum seperation value
        if maxSeriesSeperation < maxSeperation
            maxSeriesSeperation = maxSeperation;
        end
         
        % PLOT THE SEPERATION DATA ON CURRENT FIGURE
        plotHandle = plot(axesHandle,DATA.timeVector,ABSseperationTimeSeries,...
            'LineStyle','-',...
            'LineWidth',1,...
            'Color',SIM.OBJECTS(ID2).colour);
        xlabel('t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel(yAxisLabel,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylim([(-0.1*maxSeriesSeperation) (1.1*maxSeriesSeperation)]);
        grid on;

    end
end
% ADD ANNOTATION THE FOR THE VISUAL HORIZON
line([SIM.TIME.startTime SIM.TIME.endTime],...
    [SIM.OBJECTS(objectNum).detectionRange SIM.OBJECTS(objectNum).detectionRange],...
    'LineStyle','--','color','bl');
text((SIM.TIME.startTime + SIM.TIME.dt),(SIM.OBJECTS(objectNum).detectionRange - 3),...
    'Visual horizon','color','bl');   % Annotate line as visual horizon

% ADD LEGEND
legend(axesHandle,legendEntries,'Location','northeast');

% ancestor(legHandle,'axes')
% SET(axesHandle)
% handleVector = [axesHandle,legHandle]
% copyobj(handleVector,figureHandle)
% copyobj(legHandle,axesHandle)

hold off;
currentFigure = currentFigure + 1;
end
% AGENT CONTROL INPUT TRAJECTORIES
function [currentFigure,figureHandle] = get_agentControlInputs(SIM,DATA,currentFigure)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

% Declare title string for figure  
titlestr = sprintf('Agent control trajectories over %ss',num2str(DATA.timeVector(end)));

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Control Inputs');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]

plotCellA = 1; plotCellWidth = 4;                                          % The width of each figure, the start of the plot

% FOR EACH AGENT'S PERSPECTIVE
for ID1 = 1:length(DATA.objectIndex)
    % GET EVALUATION OBJECT
    evalObject = DATA.objectIndex{ID1};
    % META PROPERTIES
    simIDvector = [SIM.OBJECTS.objectID];
    simObject   = SIM.OBJECTS(simIDvector == evalObject.objectID);
    simName     = simObject.name;
    simObjectID = simObject.objectID;
    
    % CHECK IF EVALUATION OBJECT IS A WAYPOINT; OMIT IF NECESSARY
    skipCondition = simObject.type ~= OMAS_objectType.agent;
    if skipCondition 
        continue
    end
    % GET THE .DATA PROPERTY
    objectDATA = evalObject.DATA;                                          % All agents must have a .DATA property
    % OBJECT IS AN AGENT
    if isempty(objectDATA)
        warning('[OUTPUT]\tProperty .DATA is empty for agent %s',simName);
        continue
    end
    
    % GET TH AGENT-LOCAL 'DATA' STRUCTURE
    objectDATA = evalObject.DATA;
    inputTrajectories = objectDATA.inputs;
    inputCount = size(inputTrajectories,1);

    % BUILD THE SUB-PLOT
    plotCellB = ID1*plotCellWidth;                                         % The end of the plot
    plotLocation = subplot(DATA.totalAgents,plotCellWidth,[plotCellA plotCellB]);
    hold on;
    
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);                                          % Append title to first subplot
    end
    legendEntries = cell(inputCount,1);
    
    % LEGEND LABELS
    numericLabelling = 1;
    if isfield(objectDATA,'inputNames') 
        if length(objectDATA.inputNames) ~= inputCount
            warning('Incorrect number of input labels for agent %s, reverting to numeric labelling',simName);
            numericLabelling = 1;
        else
            % FETCH THE INPUT LABEL
            for entry = 1:size(inputTrajectories,1)
                legendEntries{entry} = objectDATA.inputNames{entry};
            end
            numericLabelling = 0;
        end
    end
    if numericLabelling
        % DEFAULT TO GENERIC NAMING
        for entry = 1:size(inputTrajectories,1)
            legendEntries{entry} = sprintf('Input %s',num2str(entry));
        end    
    end
    % Y AXIS LABEL
    ystring = sprintf('Control Inputs \n %s [ID:%s]',simName,num2str(simObjectID));
    
    %% PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,DATA.timeVector,inputTrajectories,...
        'LineStyle','-',...
        'LineWidth',DATA.figureProperties.LineWidth);
    xlabel('t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    ylabel(ystring,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    grid on; box on; 
    set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
          'fontWeight',DATA.figureProperties.fontWeight);
    % ADD LEGEND
    legend(legendEntries,'Location','southeast');
    % Move to next subplot location       
    plotCellA = plotCellA + plotCellWidth; 
end

hold off;
% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'inputs.fig'));
% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end

%% GENERAL OVERVIEW FIGURES
% GET THE STANDARD 3D TRAJECTORY figureProperties
function [currentFigure,figureHandle] = get_isometricFigure(SIM,DATA,currentFigure)

% FIGURE TITLE
titlestr = sprintf('System 3D trajectories over a period of %ss',num2str(DATA.timeVector(end)));

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric view');
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

legendCounter = 1; legendEntries = cell(DATA.totalObjects,1);
hold on; grid on; box on;
for ID1 = 1:DATA.totalObjects    
    % GET OBJECT OVERVIEW DATA
    legendString = sprintf('[ID:%s] %s',num2str(SIM.OBJECTS(ID1).objectID),SIM.OBJECTS(ID1).name);
    legendEntries(legendCounter) = {legendString};
    % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
    [finalStates] = OMAS_getTrajectoryData(DATA,ID1,'last');
    finalPosition = finalStates(1:3,:);
    % PLOT THE TERMINAL POSITIONS
    plot3(finalPosition(1),finalPosition(2),finalPosition(3),...
          'Marker',SIM.OBJECTS(ID1).symbol,...
          'MarkerSize',DATA.figureProperties.MarkerSize,...
          'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
          'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
          'LineWidth',DATA.figureProperties.LineWidth,...
          'LineStyle',DATA.figureProperties.LineStyle,...
          'Color',SIM.OBJECTS(ID1).colour); 
    legendCounter = legendCounter + 1;
end
legend('off')
hold on;
for ID1 = 1:DATA.totalObjects 
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    [objectStates] = OMAS_getTrajectoryData(DATA,ID1,'valid');
    positions = objectStates(1:3,:);
    plot3(positions(1,:),positions(2,:),positions(3,:),...
          'LineStyle',DATA.figureProperties.LineStyle,...
          'LineWidth',DATA.figureProperties.LineWidth,...
          'Color',SIM.OBJECTS(ID1).colour);
end
legend(legendEntries,'location','northeast');
title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);  
ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
view([45 36]);
set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
      'fontWeight',DATA.figureProperties.fontWeight);
grid on; 
hold off;
% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'isometricFigure.fig'));        % As matlab figure
saveas(figureHandle,strcat(SIM.outputPath,'isometricFigure'),'epsc');      % As postscript

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
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
for step = 1:length(DATA.timeVector)
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
            tailLength = DATA.figureProperties.tailLength;
            if activeStates < tailLength
                tailLength = activeStates;
            end
           % OBJECT HAS BEEN REMOVED
           markerPosition = objectStates(1:3,end);
           tailTrajectories = objectStates(1:3,((end+1)-tailLength):end);
        elseif step > DATA.figureProperties.tailLength
           % ONLY THE REQUESTED TAIL LENGTH
           markerPosition = objectStates(1:3,step);
           tailTrajectories = objectStates(1:3,(step - DATA.figureProperties.tailLength):step);
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

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'fourView.fig'));               % As matlab figure
saveas(figureHandle,strcat(SIM.outputPath,'fourView'),'epsc');             % As postscript

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE 3D TRAJECTORY TRAILS AS A GIF
function [currentFigure,figureHandle] = get_isometricGif(SIM,DATA,currentFigure)
% This function generates an animated .gif representing the object
% trajectories over the complete timeseries.

% OUTPUT FILE
fileName = strcat(SIM.outputPath,'isometricAnimation','.gif');

% DATA CONTAINERS
xData = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
yData = xData; zData = xData;
xMarker = NaN(1,SIM.totalObjects,SIM.TIME.numSteps); 
yMarker = xMarker; zMarker = xMarker;
% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.numSteps
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerPosition = objectStates(1:3,step);                       % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(1:3,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(1:3,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerPosition = objectStates(1:3,end);                        % Step position
            tailTrace = objectStates(1:3,((activeSteps + 1) - tailStepLength):end);
        end
        
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        xMarker(1,indexValue,step) = markerPosition(1);
        yMarker(1,indexValue,step) = markerPosition(2);
        zMarker(1,indexValue,step) = markerPosition(3); % Coordinate by object by step
        % BUILD MATICES OF TAIL COORDINATES
        traceData = tailTrace';
        xData(1:size(traceData,1),indexValue,step) = traceData(:,1); % Full X coordinate set
        yData(1:size(traceData,1),indexValue,step) = traceData(:,2); % Full X coordinate set        
        zData(1:size(traceData,1),indexValue,step) = traceData(:,3); % Full X coordinate set
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Animation');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');
set(figureHandle,'Visible','off');

% CONFIGURE VIDEO SETTINGS
fps = 50;                                                                  % Designed 50 frames per second (fps)
numFrames = SIM.TIME.simTime*fps;
stepsPerFrame = floor(SIM.TIME.numSteps/numFrames);                        % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end
F(numFrames) = struct('cdata',[],'colormap',[]);
% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(DATA.timeVector(end)));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1;
for step = 1:SIM.TIME.numSteps
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    
    % BUILD THE FRAME 
    if frame == 1
        % INITIAL TRACE PLOT
        hold on;
        traceHandle = plot3(ax,xData(:,:,step),yData(:,:,step),zData(:,:,step));
        for ind = 1:SIM.totalObjects
            % INITIAL MARKER PLOT
            markerHandle(ind) = plot3(ax,xMarker(1,ind,step),yMarker(1,ind,step),zMarker(1,ind,step));
            % ASSIGN MARKER PROPERTIES
            markerHandle(ind).Marker = SIM.OBJECTS(ind).symbol;
            markerHandle(ind).MarkerSize = DATA.figureProperties.MarkerSize;
            markerHandle(ind).MarkerFaceColor = SIM.OBJECTS(ind).colour;
            markerHandle(ind).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
            markerHandle(ind).Color = SIM.OBJECTS(ind).colour;
            % ASSIGN TRACE PROPERTIES
            traceHandle(ind).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(ind).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(ind).Color = SIM.OBJECTS(ind).colour;
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        grid on; box on; hold off;
        legend(markerHandle,DATA.figureProperties.legendEntries);
        view([-45 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            % UPDATE MARKER DATA
            set(markerHandle(entity),'XData',xMarker(:,entity,step));
            set(markerHandle(entity),'YData',yMarker(:,entity,step));
            set(markerHandle(entity),'ZData',zMarker(:,entity,step));
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity,step));
            set(traceHandle(entity),'YData', yData(:,entity,step));
            set(traceHandle(entity),'ZData', zData(:,entity,step));
        end
    end
    drawnow();
    
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    im = frame2im(F(frame));
    [imind,cm] = rgb2ind(im,256);
    
    % APPEND THE FRAMES TO GIF
    if frame == 1
        imwrite(imind,cm,fileName,'gif', 'Loopcount',inf,'DelayTime',(1/fps));
    else
        imwrite(imind,cm,fileName,'gif','WriteMode','append','DelayTime',(1/fps));
    end
    frame = frame + 1;                  % Move to next frame
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE 3D TRAJECTORY TRAIS AS A VIDEO
function [currentFigure,figureHandle] = get_isometricVideo(SIM,DATA,currentFigure)
% This function generates a video of the complete simulation trajectories
% over the complete series.

% OUTPUT FILE
fileName = strcat(SIM.outputPath,'isometricVideo','.avi');

% DATA CONTAINERS
xData = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
yData = xData; zData = xData;
xMarker = NaN(1,SIM.totalObjects,SIM.TIME.numSteps); 
yMarker = xMarker; zMarker = xMarker;
% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.numSteps
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerPosition = objectStates(1:3,step);                       % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(1:3,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(1:3,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerPosition = objectStates(1:3,end);                        % Step position
            tailTrace = objectStates(1:3,((activeSteps + 1) - tailStepLength):end);
        end
        
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        xMarker(1,indexValue,step) = markerPosition(1);
        yMarker(1,indexValue,step) = markerPosition(2);
        zMarker(1,indexValue,step) = markerPosition(3); % Coordinate by object by step
        % BUILD MATICES OF TAIL COORDINATES
        traceData = tailTrace';
        xData(1:size(traceData,1),indexValue,step) = traceData(:,1); % Full X coordinate set
        yData(1:size(traceData,1),indexValue,step) = traceData(:,2); % Full X coordinate set        
        zData(1:size(traceData,1),indexValue,step) = traceData(:,3); % Full X coordinate set
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Video Timelapse');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');
set(figureHandle,'Visible','off');

% CONFIGURE VIDEO SETTINGS
fps = 50;                                                                  % Designed 50 frames per second (fps)
numFrames = SIM.TIME.simTime*fps;
stepsPerFrame = floor(SIM.TIME.numSteps/numFrames);                        % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end
vidFile = VideoWriter(fileName);
vidFile.FrameRate = fps;                                                   % Set the videos fps to match the sample rate
open(vidFile);
% FRAME SETTINGS
frameSet = moviein(numFrames);                                             % Pre-allocate video object for frames
F(numFrames) = struct('cdata',[],'colormap',[]);
% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(DATA.timeVector(end)));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1; viewPoint = -45;
for step = 1:SIM.TIME.numSteps
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % INITIAL TRACE PLOT
        traceHandle = plot3(ax,xData(:,:,step),yData(:,:,step),zData(:,:,step));
        for ind = 1:SIM.totalObjects
            % INITIAL MARKER PLOT
            markerHandle(ind) = plot3(ax,xMarker(1,ind,step),yMarker(1,ind,step),zMarker(1,ind,step));
            % ASSIGN MARKER PROPERTIES
            markerHandle(ind).Marker = SIM.OBJECTS(ind).symbol;
            markerHandle(ind).MarkerSize = DATA.figureProperties.MarkerSize;
            markerHandle(ind).MarkerFaceColor = SIM.OBJECTS(ind).colour;
            markerHandle(ind).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
            markerHandle(ind).Color = SIM.OBJECTS(ind).colour;
            % ASSIGN TRACE PROPERTIES
            traceHandle(ind).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(ind).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(ind).Color = SIM.OBJECTS(ind).colour;
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        grid on; box on; hold off;
        legend(markerHandle,DATA.figureProperties.legendEntries);
        view([viewPoint 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            % UPDATE MARKER DATA
            set(markerHandle(entity),'XData',xMarker(:,entity,step));
            set(markerHandle(entity),'YData',yMarker(:,entity,step));
            set(markerHandle(entity),'ZData',zMarker(:,entity,step));
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity,step));
            set(traceHandle(entity),'YData', yData(:,entity,step));
            set(traceHandle(entity),'ZData', zData(:,entity,step));
        end
        view([viewPoint 50]);
        viewPoint = viewPoint + 0.05;
    end
    drawnow();
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    writeVideo(vidFile,F(frame));       % Write frame to video
    frame = frame + 1;                  % Move to next frame
end
close(vidFile);
% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end

%% TIMING FIGURES
% COMPUTATION TIME FIGURES
function [currentFigure,figureHandle] = get_computationTimes(SIM,DATA,currentFigure)
% This function gets the computation time figures for all agents using the
% DATA.objectIndex.DATA fields.
% INPUTS:
% SIM  - Local copy of the META structure
% DATA - The output data structure
% currentFigure - The current figure number
% OUTPUTS:
% currentFigure - Updated figure number
% figureHandle - Handle to the created figure

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Computation Time Series');
set(figureHandle, 'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

for ind = 1:numel(SIM.OBJECTS)
    % GET OBJECT DATA
    entity = DATA.objectIndex{ind};
    % TEST FOR THE 'DATA' FIELD
    if ~isprop(entity,'DATA') || isempty(entity.DATA)
        continue
    end
    % GET THE EQUIVALENT SIM OBJECT
    IDvector = [SIM.OBJECTS.objectID];
    SIMobject = SIM.OBJECTS(IDvector == entity.objectID);
    % GET THE TIMESERIES DATA
    dt_timeSeries = entity.DATA.algorithm_dt*1E3; 
    % LEGEND LABEL
    displayString = sprintf('[ID:%s] %s',num2str(SIMobject.objectID),SIMobject.name);
    % GENERATE FIGURES
    plot(DATA.timeVector,dt_timeSeries,...
         'Color',SIMobject.colour,...
         'LineWidth',DATA.figureProperties.LineWidth,...
         'DisplayName',displayString);
    hold on;
end
grid on; box on;
title('Agent Computation times','fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
xlabel('Simulation Time(s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
ylabel('Computation Time(ms)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
% legend('Location','northeast');
legend('Location','best');
set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
      'fontWeight',DATA.figureProperties.fontWeight);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'computationTimes.fig'));               % As matlab figure
saveas(figureHandle,strcat(SIM.outputPath,'computationTimes'),'epsc');             % As postscript
% UPDATE THE FIGURE NUMBER
if ~isempty(figureHandle)
    currentFigure = currentFigure + 1;
end
end