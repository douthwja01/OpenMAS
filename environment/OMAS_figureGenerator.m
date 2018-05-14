%% OPENMAS FIGURE GENERATOR (OMAS_figureGenerator.m) %%%%%%%%%%%%%%%%
% This script contains an index of the figures that can be requested from
% the simulation using the SIM.figures attribute.

% Author: James A. Douthwaite 10/10/2016

%% INDEX FUNCTION
function [figureNumber] = OMAS_figureGenerator(SIM,DATA,figureNumber,figureLabel)
% INPUTS:
% DATA         - The simulation output DATA structure
% figureNumber - The current figure number
% figureFlag   - The figure identifier
% OUTPUT:
% figureNumber - The updated figure number

% PARAMETERISE ALL THE FIGURES BY CURRENT CONFIGURATION
[DATA.figureProperties] = OMAS_figureProperties(SIM.OBJECTS,DATA.figureProperties);

% DETERMINE WHICH FIGURE IS TO BE GENERATED
switch upper(char(figureLabel))
    case 'ALL'
        fprintf('[%s]\tAll output figures requested.\n',SIM.phase);
        % MOVE THROUGH THE COMPLETE FIGURE VECTOR
        figureVector = {'EVENTS','COLLISIONS','TRAJECTORIES','SEPARATIONS',...
                        'CLOSEST','INPUTS','PLAN','FIG','4VIEW','GIF','AVI','TIMES'};
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
        for objectNum = 1:DATA.totalObjects
            [figureNumber,figureSet(objectNum)] = get_objectTrajectory(SIM,DATA,figureNumber,objectNum);
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = OMAS_figureTabUtility(figureSet,'OpenMAS Trajectory Overview');
        set(windowHandle,'Position',DATA.figureProperties.windowSettings);       % Maximise the figure in the tab
        savefig(windowHandle,strcat(SIM.outputPath,'globalTrajectory_overview')); % Save the output figure
    case 'SEPARATIONS' 
        fprintf('[%s]\tGenerating object trajectory separation figure.\n',SIM.phase);
        % INPUT HANDLING
        if SIM.totalObjects - SIM.totalWaypoints <= 1
            warning('There must be at least two collidable objects to plot seperation data.\n');
            figureSet = [];
            return
        else
            figureSet = [];
            for agentNum = 1:DATA.totalAgents
                [figureNumber,figureSet(agentNum)] = get_agentSeparationFigure(SIM,DATA,figureNumber,agentNum);
            end
            % ASSEMBLE TABBED FIGURE
            windowHandle = OMAS_figureTabUtility(figureSet,'OpenMAS Separation Overview');
            set(windowHandle,'Position', DATA.figureProperties.windowSettings);  % Maximise the figure in the tab
            savefig(windowHandle,strcat(SIM.outputPath,'separations_overview')); % Save the output figure
        end
    case 'CLOSEST' 
        fprintf('[%s]\tGenerating object trajectory seperation figure.\n',SIM.phase);
        [figureNumber,~] = get_minimumSeparationFigure(SIM,DATA,figureNumber);
    case 'INPUTS'
        fprintf('[%s}\tGenerating control input figure.\n',SIM.phase);
        [figureNumber,~] = get_agentControlInputs(SIM,DATA,figureNumber);
    case 'PLAN'
        fprintf('[%s}\tGenerating top-down 2D(plan) figure.\n',SIM.phase);
        [figureNumber,~] = get_topDownView(SIM,DATA,figureNumber);
    case 'FIG' 
        fprintf('[%s]\tGenerating isometric trajectory figure.\n',SIM.phase);
        [figureNumber,~] = get_isometricFigure(SIM,DATA,figureNumber);
    case '4VIEW'
        fprintf('[%s]\tGenerating four-view panel & gif file.\n',SIM.phase);
        [figureNumber,~] = get_fourViewPanel(SIM,DATA,figureNumber);
    case 'GIF'
        fprintf('[%s]\tGenerating trajectory gif file.\n',SIM.phase);
        [figureNumber,~] = get_isometricGif(SIM,DATA,figureNumber);
    case 'AVI'
        fprintf('[%s]\tGenerating trajectory avi file.\n',SIM.phase);
        [figureNumber,~] = get_isometricAvi(SIM,DATA,figureNumber);
    case 'TIMES'
        fprintf('[%s]\tGenerating computation timeseries figure.\n',SIM.phase);
        [figureNumber,~] = get_computationTimes(SIM,DATA,figureNumber);
    case 'NONE'
        fprintf('[%s]\tFigure output supressed.\n',SIM.phase);
    otherwise
        warningStr = sprintf('[ERROR] Did not recognise figure request: "%s"',char(figureLabel));
        warning(warningStr);
end
end


%% FIGURE GENERATION FUNCTIONS 
% ////////////////// COLLISION DATA FIGURE GENERATION /////////////////////
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
figurePath = strcat(SIM.outputPath,'eventOverview');
figureHandle = figure('Name','OpenMAS Event Overview');
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
xlabel('t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
ylabel('Number of Events','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
set(h,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
set(h,'Color',DATA.figureProperties.axesColor);
xlim([-SIM.TIME.dt,SIM.TIME.endTime])
legend(h,eventFields,'Location','northeastoutside')
grid on;

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
tbox = annotation(figureHandle,'textbox',[columnXOffset textboxYOffset columnWidth textBoxHeight],'Units','Normalized','String',overviewInfo); 
set(tbox,'BackgroundColor','w',...
               'EdgeColor','k',...
                  'Margin', 5,...
                'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
            
% BUILD THE NEIGHBOURING TEXT BLOCKS
h = subplot(subHandles(3));   
axis off;  
% CREATE THE COLLISION DATA CONTENTS
agentSummary    = sprintf('Agent Collisions: \n %s/%s (%3.0f%%)',num2str(DATA.collisions),num2str(DATA.totalAgents),DATA.collisionPercentage);
waypointSummary = sprintf('Waypoints: \n %s/%s (%3.0f%%)',num2str(DATA.waypointsAchieved),num2str(DATA.totalWaypoints),DATA.waypointPercentage);
statsInfo = {'Scenario Summary','',agentSummary,waypointSummary};


% CREATE WAYPOINT DATA OVERVIEW
textboxYOffset = 0.11;
textBoxHeight = 0.3;
% tbox = annotation(figureHandle,'textbox',[0.65 0.11 .32 .4],'Units','Normalized','String',statsInfo); 
tbox = annotation(figureHandle,'textbox',[columnXOffset textboxYOffset columnWidth textBoxHeight],'Units','Normalized','String',statsInfo); 
set(tbox,'BackgroundColor','w',...
               'EdgeColor','k',...
                  'Margin', 5,...
                'FontSize',DATA.figureProperties.axisFontSize,...
              'fontWeight',DATA.figureProperties.fontWeight);
          
title(subplot(subHandles(1)),sprintf('Event occurance over a period of %ss',num2str(SIM.TIME.endTime)),...
       'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize,'FontSmoothing','on');
% suptitle(sprintf('Event occurance over the %ss period)',num2str(DATA.timeVector(end))));

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

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
    figureHandle = figure('Name',tabStr);
    set(figureHandle,'Position', DATA.figureProperties.windowSettings);       % [x y width height]
    set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
    
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
    radiusA = SIM.OBJECTS(collisionEvent.objectID_A).radius;
    [X,Y,Z] = sphere(40);
    Xa = X.*radiusA + collisionEvent.state_A(1);
    Ya = Y.*radiusA + collisionEvent.state_A(2);
    Za = Z.*radiusA + collisionEvent.state_A(3);    
    mesh(Xa,Ya,Za,...
        'FaceColor',SIM.OBJECTS(collisionEvent.objectID_A).colour,...
        'FaceAlpha',0.1,...
        'LineWidth',DATA.figureProperties.LineWidth,...
        'EdgeAlpha',0.2,...
        'edgecolor',DATA.figureProperties.MarkerEdgeColor);  % Obstacle

    % Build object B volume representation
    radiusB = SIM.OBJECTS(collisionEvent.objectID_B).radius;
    Xb = X.*radiusB + collisionEvent.state_B(1);
    Yb = Y.*radiusB + collisionEvent.state_B(2);
    Zb = Z.*radiusB + collisionEvent.state_B(3);
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
    xlabel('x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    ylabel('y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    zlabel('z(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
          'fontWeight',DATA.figureProperties.fontWeight);
    set(gca,'Color',DATA.figureProperties.axesColor);
    hold off;
    
    figureSet = vertcat(figureSet,figureHandle);
end

% ASSEMBLE TABBED FIGURE
tabbedFigureHandle = OMAS_figureTabUtility(figureSet,'OpenMAS Collision Overview');

% SAVE THE OUTPUT FIGURE
savefig(tabbedFigureHandle,strcat(SIM.outputPath,'collisionOverview.fig'));

currentFigure = currentFigure + 1;
end

% ////////////////// TRAJECTORY DATA FIGURE GENERATION ////////////////////
% GET THE GLOBAL TRAJECTORY DATA FOR AN INDIVIDUAL OBJECT
function [currentFigure,figureHandle] = get_objectTrajectory(SIM,DATA,currentFigure,objectNum)

% Declare title string for figure  
titlestr = sprintf('Global trajectory data for %s over a period of %ss',SIM.OBJECTS(objectNum).name,num2str(SIM.TIME.endTime));

% CONFIGURE THE PLOT ATTRIBUTES
objectLabel = sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(objectNum).objectID),SIM.OBJECTS(objectNum).name);
figureHandle = figure('Name',objectLabel);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'Visible','off');
plotCellWidth = 4; plotCellA = 1;                                          % The width of each figure, the start of the plot

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[objectStates] = OMAS_getTrajectoryData(DATA,objectNum);              % Get the object data
% STATE NAME VECTOR
stateTags = {'x(m)','y(m)','z(m)',...
             'u(m/s)','v(m/s)','w(m/s)',...
             'q0','q1','q2','q3'};

setappdata(gcf, 'SubplotDefaultAxesLocation', [0.1,0.1,0.85,0.82]);         

% FOR EACH OBJECTS PERSPECTIVE
stateVectorLength = size(objectStates,1);
for stateNum = 1:stateVectorLength
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    stateTrajectory = objectStates(stateNum,:);
        
    % CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
    plotCellB = stateNum*plotCellWidth;                                    % The end of the plot
    plotLocation = subplot(stateVectorLength,plotCellWidth,[plotCellA plotCellB]);
    
    % PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,DATA.timeVector(1:SIM.TIME.endStep),stateTrajectory(:,1:SIM.TIME.endStep),...
         'LineStyle','-',...
         'LineWidth',DATA.figureProperties.LineWidth,...
         'Color','b');
    ylabel(plotLocation,stateTags{stateNum},'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
    set(plotLocation,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
    set(plotLocation,'Color',DATA.figureProperties.axesColor);
    grid on; box on; grid minor;
    
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if stateNum == 1
        title(plotLocation,titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize,'FontSmoothing','on');                                                   % Append title to first subplot
    end
    % Prevent overlap of x-label
    if stateNum == stateVectorLength
        xlabel(plotLocation,'t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
    else
        set(plotLocation,'XTickLabel',[]); 
    end
    % Move to next subplot location 
    plotCellA = plotCellA + plotCellWidth;
end
hold off; 

% SAVE THE OUTPUT FIGURE
fileName = strcat(SIM.outputPath,sprintf('globalTrajectory_[ID%s] %s',num2str(SIM.OBJECTS(objectNum).objectID),SIM.OBJECTS(objectNum).name));
savefig(figureHandle,fileName);                                            % As matlab figure 

% ITERATE PLOT
currentFigure = currentFigure + 1;
end
% AGENT CONTROL INPUT TRAJECTORIES
function [currentFigure,figureHandle] = get_agentControlInputs(SIM,DATA,currentFigure)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

% Declare title string for figure  
titlestr = sprintf('Agent control trajectories over a period of %ss',num2str(SIM.TIME.endTime));
figurePath = strcat(SIM.outputPath,'inputs');

% FIGURE META PROPERTIES
figureHandle = figure('Name','OpenMAS Control Inputs');
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.88, 0.85]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
plotCellA = 1; plotCellWidth = 4;                                          % The width of each figure, the start of the plot

% LOGICALLY GET THE AGENT'S FROM THE OBJECT INDEX
agentObjectIndex = DATA.objectIndex([SIM.OBJECTS.type] == OMAS_objectType.agent); % The agents themselves

% FOR EACH AGENT'S PERSPECTIVE
for ID1 = 1:SIM.totalAgents
    % GET EVALUATION OBJECT
    evalAgent = agentObjectIndex{ID1};
    % OBJECT IS AN AGENT
    if isempty(evalAgent.DATA)
        warning('[OUTPUT]\tProperty .DATA is empty for agent %s',evalAgent.name);
        continue
    end
    
    % GET TH AGENT-LOCAL 'DATA' STRUCTURE
    objectDATA = evalAgent.DATA;
    inputTrajectories = objectDATA.inputs;
    inputCount = size(inputTrajectories,1);

    % BUILD THE SUB-PLOT
    plotCellB = ID1*plotCellWidth;                                         % The end of the plot
    plotLocation = subplot(DATA.totalAgents,plotCellWidth,[plotCellA plotCellB]);
    hold on;
    % LEGEND LABELS
    legendEntries = cell(inputCount,1);   
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
    ystring = sprintf('Inputs \n [ID:%s]',num2str(evalAgent.objectID));%,evalAgent.name);
    
    % PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,DATA.timeVector(:,1:SIM.TIME.endStep),inputTrajectories(:,1:SIM.TIME.endStep),...
        'LineStyle','-','LineWidth',DATA.figureProperties.LineWidth);
    
    ylabel(plotLocation,ystring,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','On');
    grid on; box on; 
    set(plotLocation,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','On');
    set(plotLocation,'Color',DATA.figureProperties.axesColor);
    
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        title(plotLocation,titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);                                          % Append title to first subplot
    end
    if ID1 == DATA.totalAgents
        xlabel(plotLocation,'t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','On');                                        % Append title to first subplot
    end
    
    % ADD LEGEND
    legend(legendEntries,'Location','eastoutside');
    % Move to next subplot location       
    plotCellA = plotCellA + plotCellWidth; 
end
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end
% GET THE AGENT SEPERATION SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = get_agentSeparationFigure(SIM,DATA,currentFigure,agentNum)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

% INPUT HANDLING
if SIM.totalObjects - SIM.totalWaypoints <= 1
    warning('There must be at least two collidable objects to plot seperation data.\n');
    figureHandle = [];
    return
end

% GENERATE THE SEPERATION DATA FOR A GIVEN AGENT
agentMETA = SIM.OBJECTS([SIM.OBJECTS.type] == OMAS_objectType.agent);      % The agent set
agentSubject = agentMETA(agentNum);                                        % The agent requested

globalIndex =  find((SIM.globalIDvector == agentSubject.objectID), 1, 'first'); % The agents position in the global series
agentLabel = sprintf('[ID-%s] %s',num2str(agentSubject.objectID),agentSubject.name);
figurePath = strcat(SIM.outputPath,sprintf('separations_%s ',agentLabel));

% FIGURE META DATA
figureHandle = figure('Name',agentLabel);                                   % Tab label
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'Visible','off');

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[agentStates] = OMAS_getTrajectoryData(DATA,globalIndex);                  % Trajectory data per objectID

% CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
legendEntries = cell(1,(DATA.totalAgents-1));
legendCounter = 1; maxSeriesSeperation = 1;                                % Initial y-axis limit

ax = axes(figureHandle);
hold on;

% MOVE THROUGH OTHER OBJECTS
for ID2 = 1:DATA.totalObjects
    % DEFINE CONDITIONS FOR PLOTTING
    plotCondition = agentSubject.objectID ~= SIM.OBJECTS(ID2).objectID && ~strcmpi(SIM.OBJECTS(ID2).type,'WAYPOINT');
    if plotCondition
        % Generate legend entry
        legendEntries(legendCounter) = {sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(ID2).objectID),num2str(SIM.OBJECTS(ID2).name))};
        % Container for the timeseries data
        ABSseperationTimeSeries = zeros(1,size(agentStates,2));
        % Get comparative state data
        [objectStates] = OMAS_getTrajectoryData(DATA,ID2);
        
        % Reset collision instance to 0
        collisionInstance = 0;
        for simStep = 1:size(ABSseperationTimeSeries,2)
            % CALCULATE THE ABSOLUTE OBSTACLE SEPERATION TIMESERIES
            collisionCondition = agentSubject.radius + SIM.OBJECTS(ID2).radius;                      % Define the collision condition (physical seperation)
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
            'LineWidth',DATA.figureProperties.LineWidth,...
            'Color',SIM.OBJECTS(ID2).colour);
        legendCounter = legendCounter + 1;                             % Increment legend entry
    end
end
% THE COLLISION CONDITION
refHandle = refline(ax,0,collisionCondition);                                % Adds a reference line with slope m and intercept b to the current axes.
set(refHandle,'color','k',...
    'LineStyle','--',...
    'LineWidth',DATA.figureProperties.LineWidth);
legendEntries{legendCounter} = 'Collision Boundary';

% ADD FINAL PLOT ATTRIBUTES
title(ax,sprintf('Object seperations for agent %s',agentLabel),'fontweight',DATA.figureProperties.fontWeight,...
                 'fontSize',DATA.figureProperties.titleFontSize);
xlabel(ax,'t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
ylabel(ax,'Seperation(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
ylim(ax,[(-0.1*maxSeriesSeperation) (1.1*maxSeriesSeperation)]);
set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
set(ax,'Color',DATA.figureProperties.axesColor);
legend(ax,legendEntries,'Location','northeastoutside');
grid on; box on;
drawnow;
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;      
end
% GET THE CLOSEST PROXIMITY FIGURES
function [currentFigure,figureHandle] = get_minimumSeparationFigure(SIM,DATA,currentFigure)
% This function generates a figure with the smallest separations observed
% in the simulation and plots them on a inter-agent basis. This function is
% irrespective of agents, or obstacles, but neglects waypoints.

% INPUT HANDLING
if SIM.totalObjects - SIM.totalWaypoints <= 1
    warning('There must be at least two collidable objects to plot seperation data.');
    figureHandle = [];
    return
end

% BUILD OBJECT SET FOR ALL NON-WAYPOINTS
objectMETA = SIM.OBJECTS([SIM.OBJECTS.type] ~= OMAS_objectType.waypoint);
[LIA,LOCB] = ismember(SIM.globalIDvector,[objectMETA.objectID]);
collidableIDs = LOCB(LIA);

% OUTPUT CONTAINER
separationTimeSeries = inf(numel(collidableIDs),numel(collidableIDs),SIM.TIME.numSteps);
for IDnumA = 1:numel(collidableIDs)
    % THE AGENT STATES
    [objectStatesA] = OMAS_getTrajectoryData(DATA,collidableIDs(IDnumA));
    
    for IDnumB = 1:numel(collidableIDs)
        if IDnumA == IDnumB 
            % OMIT SEPERATIONS BETWEEN ITSELF
            continue
        else
            % GET THE AGENT STATE TIMESERIES
            objectStatesB = OMAS_getTrajectoryData(DATA,collidableIDs(IDnumB));
            centroidSeparations = objectStatesB(1:3,:) - objectStatesA(1:3,:);  % seperation of the centroids
            centroidSeparations = sqrt(sum(centroidSeparations.^2,1));
            % STORE IN SEPERATION TIMESERIES MATRIX
            separationTimeSeries(IDnumB,IDnumA,:) = centroidSeparations;
        end
    end
end

% GET THE MINIMUM SEPERATIONS FOR EACH AGENT
% If the number of obstacle/agents exceeds 10 objects, we are only
% interested in the 10 interactions that came the closest to collision.
maxDisplayObjects = inf;

% MIN AND MAXIMUM SEPERATION MATRICES [IDA by IDB]
minABAxes = collidableIDs;
minABMatrix = min(separationTimeSeries,[],3);                          % Minimum separations over the timeseries
maxABMatrix = max(separationTimeSeries,[],3);

maxABMatrix(isinf(maxABMatrix)) = NaN;

[h_maxVals,~] = max(maxABMatrix,[],2);       % The minimum seperations, horezontal indicies
maxProximity = max(h_maxVals);               % Maximal seperation

% MIN-MATRIX > MINIMUM SEPERATIONS [IDA by IDB]


% minABMatrix is a matrix of closest interactions between all objects A
% vs B.
[h_val,h_ind] = min(minABMatrix,[],2); % The minimum seperations, horezontal indicies

% DETERMINE UNIQUE INTERACTIONS
% [ord,ord_ind] = sort(h_val);           % Sort minimums by magnitude
% [~,unique_ind,~] = unique(ord,'rows');

% uniqueIndOrder = ord_ind(unique_ind);  % Unique, ordered proximities

% % THE ID's OF THE CLOSEST (UNIQUE) OBJECTS
closestIDpairs = [minABAxes',h_ind];
% orderedClosestPairs = closestIDpairs(uniqueIndOrder,:);
% % GET THE LIMITED NUMBER OF ID PAIRS TO PLOT
% if size(orderedClosestPairs,1) > maxDisplayObjects
%     plotableInteractions = orderedClosestPairs(1:maxDisplayObjects,:);
% else
%     plotableInteractions = orderedClosestPairs;
% end

plotableInteractions = closestIDpairs; % <<< PLOT ALL, CLOSEST INTERACTIONS, BETWEEN OBJECTS

% FIGURE META PROPERTIES
figurePath = strcat(SIM.outputPath,'minimumSeparations');   
figureHandle = figure('Name','OpenMAS point of closest approach'); 
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 

% BEGIN GENERATING THE FIGURE
ax = axes(figureHandle);
hold on;
legendEntries = {}; legendCounter = 1;
% FOR EACH INTERACTION GET THE SEPERATION DATA
for interaction = 1:size(plotableInteractions,1)
    % GET THE INTERACTION IDs
    ID1 = plotableInteractions(interaction,1);
    ID2 = plotableInteractions(interaction,2);
    
    % PLOT A LINE WITH THE COMBINED RADIUS
    criticalLimit = SIM.OBJECTS(ID1).radius + SIM.OBJECTS(ID2).radius;
    
    % BUILD LEGEND STRING
    label1 = sprintf('[ID:%s]%s',num2str(SIM.OBJECTS(ID1).objectID),num2str(SIM.OBJECTS(ID1).name));
    label2 = sprintf('[ID:%s]%s',num2str(SIM.OBJECTS(ID2).objectID),num2str(SIM.OBJECTS(ID2).name));
    
    legendEntries{legendCounter} = strcat(label1,' - ',label2);
%     legendEntries = horzcat(legendEntries,label1,label2);

    % GET THE SEPERATION TIMESERIES TO PLOT
    sepTimeSeries = squeeze(separationTimeSeries(ID1,ID2,:));
    % PLOT THE SEPERATIONS
    lineHandles(interaction) = plot(ax,DATA.timeVector,sepTimeSeries');
    set(lineHandles(interaction),...
        'color',SIM.OBJECTS(ID1).colour,...
        'LineStyle','-',...
        'LineWidth',DATA.figureProperties.LineWidth);
    % INCREMENT THE LEGEND COUNTER
    legendCounter = legendCounter + 1;                             % Increment legend entry
end

% THE COLLISION CONDITION
refHandle = refline(ax,0,criticalLimit);                                % Adds a reference line with slope m and intercept b to the current axes.
set(refHandle,'color','k','LineStyle','--','LineWidth',DATA.figureProperties.LineWidth);


% ADD FINAL PLOT ATTRIBUTES
titleString = sprintf('Global minimum separations over a period of %ss',num2str(SIM.TIME.endTime));
title(ax,titleString,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.titleFontSize,'fontSmoothing','On');
xlabel(ax,'t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'fontSmoothing','On');
ylabel(ax,'Separation(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'fontSmoothing','On');
set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'fontSmoothing','On');
set(ax,'Color',DATA.figureProperties.axesColor);
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
ylim(ax,[0 (0.5*maxProximity)]);
grid on; box on;
drawnow;
hold off;

% legendEntries{legendCounter} = 'Collision Boundary';
% legHandle = legend(ax,legendEntries,'Location','northeastoutside');

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;    
end

% /////////////////////// GENERAL OVERVIEW FIGURES ////////////////////////
% GET THE STANDARD 2D TRAJECTORY DATA.figureProperties [ UPDATED ]
function [currentFigure,figureHandle] = get_topDownView(SIM,DATA,currentFigure)

% FIGURE TITLE
titlestr = sprintf('Object trajectories over a period of %ss',num2str(SIM.TIME.endTime));
figurePath = strcat(SIM.outputPath,'2DplanFigure');

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS aerial view');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW

legendCounter = 1; legendEntries = cell(DATA.totalObjects,1);
ax = gca;
hold on; grid on; box on;
for ID1 = 1:DATA.totalObjects    
    % GET OBJECT OVERVIEW DATA
    legendString = sprintf('[ID:%s] %s',num2str(SIM.OBJECTS(ID1).objectID),SIM.OBJECTS(ID1).name);
    legendEntries(legendCounter) = {legendString};
    % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
    [finalStates] = OMAS_getTrajectoryData(DATA,ID1,'last');
    finalPosition = finalStates(1:3,:);
    % DISPLAY THE OBJECT
    if isstruct(SIM.OBJECTS(ID1).patch)
        patch(ax,'Vertices',SIM.OBJECTS(ID1).patch.vertices*SIM.OBJECTS(ID1).R_GB + finalPosition',...
            'Faces',SIM.OBJECTS(ID1).patch.faces,...
            'FaceColor',SIM.OBJECTS(ID1).colour,...
            'EdgeColor',DATA.figureProperties.EdgeColor,...
            'EdgeAlpha',DATA.figureProperties.EdgeAlpha,...  
            'FaceLighting',DATA.figureProperties.FaceLighting,...
            'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties
    else
        % PLOT THE TERMINAL POSITIONS
        plot3(ax,finalPosition(1),finalPosition(2),finalPosition(3),...
              'Marker',SIM.OBJECTS(ID1).symbol,...
              'MarkerSize',DATA.figureProperties.MarkerSize,...
              'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
              'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
              'LineWidth',DATA.figureProperties.LineWidth,...
              'LineStyle',DATA.figureProperties.LineStyle,...
              'Color',SIM.OBJECTS(ID1).colour); 
    end  
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
      
    % ADD THE INTIAL POSITION REFERENCE POSITION 
    if SIM.OBJECTS(ID1).type == OMAS_objectType.agent
        initialPlotPosition = positions(1:2,1) - SIM.OBJECTS(ID1).radius; % The first position of the object
        rectangle('Position',[initialPlotPosition(1) initialPlotPosition(2) (2*SIM.OBJECTS(ID1).radius) (2*SIM.OBJECTS(ID1).radius)],...
                  'Curvature',[1 1],...
                  'FaceColor',SIM.OBJECTS(ID1).colour,...
                  'EdgeColor',DATA.figureProperties.EdgeColor,...
                  'LineWidth',DATA.figureProperties.PatchLineWidth)  
    end  
end
legend(legendEntries,'location','northeastoutside');
% legend(ax,legendEntries,'location','best');
title(ax,titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
xlabel(ax,'x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);  
ylabel(ax,'y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
xlim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
ylim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
zlim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');% axis equal
set(ax,'Color',DATA.figureProperties.axesColor);
axis square;     
view([0 90]);
% set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE STANDARD 3D TRAJECTORY DATA.figureProperties [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricFigure(SIM,DATA,currentFigure)

% FIGURE TITLE
titlestr = sprintf('Object trajectories over a period of %ss',num2str(SIM.TIME.endTime));
figurePath = strcat(SIM.outputPath,'isometricFigure');

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric view');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW

legendCounter = 1; legendEntries = cell(DATA.totalObjects,1);
ax = gca;
hold on; grid on; box on;
for ID1 = 1:DATA.totalObjects    
    % GET OBJECT OVERVIEW DATA
    legendString = sprintf('[ID:%s] %s',num2str(SIM.OBJECTS(ID1).objectID),SIM.OBJECTS(ID1).name);
    legendEntries(legendCounter) = {legendString};
    % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
    [finalStates] = OMAS_getTrajectoryData(DATA,ID1,'last');
    finalPosition = finalStates(1:3,:);
    % DISPLAY THE OBJECT
    if isstruct(SIM.OBJECTS(ID1).patch)
        patch(ax,'Vertices',SIM.OBJECTS(ID1).patch.vertices*SIM.OBJECTS(ID1).R_GB + finalPosition',...
            'Faces',SIM.OBJECTS(ID1).patch.faces,...
            'FaceColor',SIM.OBJECTS(ID1).colour,...
            'EdgeColor',DATA.figureProperties.EdgeColor,...
            'EdgeAlpha',DATA.figureProperties.EdgeAlpha,...  
            'FaceLighting',DATA.figureProperties.FaceLighting,...
            'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties            % Patch properties
    else
        % PLOT THE TERMINAL POSITIONS
        plot3(ax,finalPosition(1),finalPosition(2),finalPosition(3),...
              'Marker',SIM.OBJECTS(ID1).symbol,...
              'MarkerSize',DATA.figureProperties.MarkerSize,...
              'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
              'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
              'LineWidth',DATA.figureProperties.LineWidth,...
              'LineStyle',DATA.figureProperties.LineStyle,...
              'Color',SIM.OBJECTS(ID1).colour); 
    end  
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
legend(legendEntries,'location','northeastoutside');
% legend(ax,legendEntries,'location','best');
title(ax,titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
xlabel(ax,'x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);  
ylabel(ax,'y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
zlabel(ax,'z(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize); 
xlim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
ylim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
zlim(ax,[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');% axis equal
set(ax,'Color',DATA.figureProperties.axesColor);
axis vis3d;     
view([-24 36]);
set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
grid on; 
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% print(figureHandle,figurePath,'-dpdf','-r0');

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE 3 PLAN + ISOMETRIC ANIMATED TRAILS FIGURE
function [currentFigure,figureHandle] = get_fourViewPanel(SIM,DATA,currentFigure)
% This function generates a four-view animated panel showing the motion of
% all objects in the simulation.

% OUTPUT FILE
fileName = strcat(SIM.outputPath,'fourView','.gif');

% DATA CONTAINERS
% We need to build a matrix of states*tailLength*IDs,step
globalTraces = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),10,DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
globalMarkers = NaN(10,SIM.totalObjects,SIM.TIME.numSteps); 

% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.endStep
    for entity = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,entity,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerStates = objectStates(:,step);                            % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(:,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(:,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerStates = objectStates(:,end);                             % Step position
            tailOrigin = ((activeSteps + 1) - tailStepLength);
            if tailOrigin > 0
                tailTrace = objectStates(:,((activeSteps + 1) - tailStepLength):end);
            else
                tailTrace = objectStates(:,(activeSteps + 1):end);
            end
        end
                
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        globalMarkers(:,entity,step) = markerStates;      
        
        % BUILD MATICES OF TAIL COORDINATES
        % We need to build a matrix of dimensions
        % [tailLength*states*IDs,step], the resulting plot matrices must be
        % of dimensions [state(1)*tailLength,IDs]
        globalTraces(1:size(tailTrace,2),:,entity,step) = tailTrace';
    end
end

% BUILD FIGURE
% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Four Viewpoint Panel');
% META FIGURE PROPERTIES
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.1, 0.90, 0.88]);
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');

% CONFIGURE VIDEO SETTINGS
fps = 50; % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.endTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end
% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(SIM.TIME.endStep/numFrames);                         % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end

F(numFrames) = struct('cdata',[],'colormap',[]);
% AXES SETTINGS
titleString = sprintf('Object planar trajectories over a period of %ss',num2str(SIM.TIME.endTime));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1;

% GET THE LAST INSTANCE OF A VALID STATE
for step = 1:SIM.TIME.endStep
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % STEP DATA
    xData = squeeze(globalTraces(:,1,:,step));
    yData = squeeze(globalTraces(:,2,:,step));
    zData = squeeze(globalTraces(:,3,:,step));
    markerStates = globalMarkers(:,:,step);
    
    % BUILD THE FRAME 
    if frame == 1
        % GENERATE THE SUB PLOTS (SPECIFIC ORIENTATIONS)
        subHandles(1) = subplot(4,4,[1 6]);%,'align'
        subHandles(2) = subplot(4,4,[3 8]);%,'align'
        subHandles(3) = subplot(4,4,[9 14]);%,'align'
        subHandles(4) = subplot(4,4,[11 16]);%'align');           
        hold on;
        for entity = 1:SIM.totalObjects
            % GET THE GLOBAL POSE AT STEP
            if isstruct(SIM.OBJECTS(entity).patch)
                % CALCULATE PATCH 'POSE'
                [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(markerStates(7:end,entity));
                objectVertices = SIM.OBJECTS(entity).patch.vertices*R2 + markerStates(1:3,entity)';
                % GENERATE MARKER HANDLES
                for handleID = 1:numel(subHandles)
                    h = subplot(subHandles(handleID));
                    % GENERATE PATCH PLOTS
                    markerHandles(handleID,entity) = patch(h,'Vertices',objectVertices,...
                        'Faces',SIM.OBJECTS(entity).patch.faces,...
                        'FaceColor',SIM.OBJECTS(entity).colour,...
                        'EdgeColor',DATA.figureProperties.EdgeColor,...
                        'EdgeAlpha',DATA.figureProperties.EdgeAlpha,... 
                        'FaceLighting',DATA.figureProperties.FaceLighting,...
                        'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties
                    hold on;
                    % ASSIGN TRACE PROPERTIES
                    traceHandles(handleID,entity) = plot3(h,xData(:,entity),yData(:,entity),zData(:,entity),...
                    'LineStyle',DATA.figureProperties.LineStyle,...
                    'LineWidth',DATA.figureProperties.LineWidth,...
                    'Color',SIM.OBJECTS(entity).colour);
                    xlabel(h,'x(m)'); ylabel(h,'y(m)'); zlabel(h,'z(m)');
                end
            else
                % GENERATE MARKER HANDLES
                for handleID = 1:numel(subHandles)
                    % OTHERWISE USE THEIR ALLOCATED SYMBOL
                    h = subplot(subHandles(handleID));
                    markerHandles(handleID,entity) = plot3(h,markerStates(1,entity),markerStates(2,entity),markerStates(3,entity),...
                    'Marker',SIM.OBJECTS(entity).symbol,...
                    'MarkerSize',DATA.figureProperties.MarkerSize,...
                    'MarkerFaceColor',SIM.OBJECTS(entity).colour,...
                    'MarkerEdgeColor',DATA.figureProperties.MarkerEdgeColor,...
                    'Color',SIM.OBJECTS(entity).colour);
                    hold on;
                    % ASSIGN TRACE PROPERTIES
                    traceHandles(handleID,entity) = plot3(h,xData(:,entity),yData(:,entity),zData(:,entity),...
                    'LineStyle',DATA.figureProperties.LineStyle,...
                    'LineWidth',DATA.figureProperties.LineWidth,...
                    'Color',SIM.OBJECTS(entity).colour);
                    xlabel(h,'x(m)'); ylabel(h,'y(m)'); zlabel(h,'z(m)');
                end
            end 
        end
        % CONFIGURE THE PLOT SET 
        set(subHandles(:),'Box','on','XGrid','on','YGrid','on','ZGrid','on');
        set(subHandles(:),'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
        set(subHandles(:),'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(subHandles(:),'Color',DATA.figureProperties.axesColor);
        set(subHandles(:),'XLimMode','manual','YLimMode','manual','ZLimMode','manual');
        set(subHandles(:),'XLim',[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        set(subHandles(:),'YLim',[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        set(subHandles(:),'ZLim',[DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        % ADD LEGEND
        if SIM.totalObjects <= 10
            %legend(subHandles(1),DATA.figureProperties.legendEntries,'Location','NorthEastOutside');
        end
        % CORRECT ORIENTATIONS
        set(subHandles(1),'View',[45 45]);
        set(subHandles(2),'View',[0 90]);
        set(subHandles(3),'View',[90 0]);
        set(subHandles(4),'View',[0 0]);
        % META TITLE
%         suptitle(titleString)
    else
        % UPDATE THE EXISTING PLOT HANDLES FOR EACH AGENT
        for entity = 1:SIM.totalObjects
            for handleID = 1:numel(subHandles)
                
                if isstruct(SIM.OBJECTS(entity).patch)
                    % GET THE GLOBAL POSE AT STEP
                    [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(markerStates(7:end,entity));
                    % BUILD MARKER
                    set(markerHandles(handleID,entity),'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + markerStates(1:3,entity)');
                else
                    % UPDATE MARKER DATA
                    set(markerHandles(handleID,entity),'XData',markerStates(1,entity));
                    set(markerHandles(handleID,entity),'YData',markerStates(2,entity));
                    set(markerHandles(handleID,entity),'ZData',markerStates(3,entity));
                end
                % UPDATE TRACE DATA
                set(traceHandles(handleID,entity),'XData', xData(:,entity));
                set(traceHandles(handleID,entity),'YData', yData(:,entity));
                set(traceHandles(handleID,entity),'ZData', zData(:,entity));
            end
        end
    end
    drawnow;
    
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
% GET THE 3D TRAJECTORY TRAILS AS A GIF  [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricGif(SIM,DATA,currentFigure)
% This function generates an animated .gif representing the object
% trajectories over the complete timeseries.

% OUTPUT FILE
filePath = strcat(SIM.outputPath,'isometricFigure.gif');

% DATA CONTAINERS
% We need to build a matrix of states*tailLength*IDs,step
globalTraces = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),10,DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
globalMarker = NaN(10,SIM.totalObjects,SIM.TIME.numSteps); 

% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.endStep
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerState = objectStates(:,step);                            % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(:,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(:,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerState = objectStates(:,end);                             % Step position
            tailTrace = objectStates(:,((activeSteps + 1) - tailStepLength):end);
        end
                
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        globalMarker(:,indexValue,step) = markerState;      
        
        % BUILD MATICES OF TAIL COORDINATES
        % We need to build a matrix of dimensions
        % [tailLength*states*IDs,step], the resulting plot matrices must be
        % of dimensions [state(1)*tailLength,IDs]
        globalTraces(1:size(tailTrace,2),:,indexValue,step) = tailTrace';
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Timelapse');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');
% set(figureHandle,'Visible','off');


% CONFIGURE VIDEO SETTINGS
fps = 50; % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.endTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end
% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(SIM.TIME.endStep/numFrames);                         % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end

F(numFrames) = struct('cdata',[],'colormap',[]);
% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1;
for step = 1:SIM.TIME.endStep
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % STEP DATA
    xData = squeeze(globalTraces(:,1,:,step));
    yData = squeeze(globalTraces(:,2,:,step));
    zData = squeeze(globalTraces(:,3,:,step));
    
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % INITIAL TRACE PLOT
        traceHandle = plot3(ax,xData,yData,zData);
        % UPDATE OBJECT HANDLES
        for entity = 1:SIM.totalObjects
            if isstruct(SIM.OBJECTS(entity).patch)
               % GET THE GLOBAL POSE AT STEP
               [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,entity,step));
               % BUILD MARKER
               markerHandle(entity) = patch(ax,'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + globalMarker(1:3,entity,step)',...
                   'Faces',SIM.OBJECTS(entity).patch.faces,...
                   'FaceColor',SIM.OBJECTS(entity).colour,...
                   'EdgeColor',DATA.figureProperties.EdgeColor,...
                   'EdgeAlpha',DATA.figureProperties.EdgeAlpha,...          
                   'FaceLighting',DATA.figureProperties.FaceLighting,...
                   'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(entity) = plot3(ax,globalMarker(1,entity,step),globalMarker(2,entity,step),globalMarker(3,entity,step));
                % ASSIGN MARKER PROPERTIES
                markerHandle(entity).Marker = SIM.OBJECTS(entity).symbol;
                markerHandle(entity).MarkerSize = DATA.figureProperties.MarkerSize;
                markerHandle(entity).MarkerFaceColor = SIM.OBJECTS(entity).colour;
                markerHandle(entity).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
                markerHandle(entity).Color = SIM.OBJECTS(entity).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(entity).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(entity).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(entity).Color = SIM.OBJECTS(entity).colour;
            % ANNOTATION WITH THE CURRENT TIME
            clockHandle = annotation('textbox',[0.02 0.05 0.15 0.04],'String','Time:','FitBoxToText','off');
            set(clockHandle,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        zlabel('z(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
        set(ax,'Color',DATA.figureProperties.axesColor);
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        set(ax,'outerposition',[0.02 0.1 1 0.88]);
        grid on; box on; hold off;
        legend(markerHandle,DATA.figureProperties.legendEntries,'Location','northeastoutside');
        view([-45 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            if isstruct(SIM.OBJECTS(entity).patch)
                % GET THE GLOBAL POSE AT STEP
                [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,entity,step));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + globalMarker(1:3,entity,step)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',globalMarker(1,entity,step));
                set(markerHandle(entity),'YData',globalMarker(2,entity,step));
                set(markerHandle(entity),'ZData',globalMarker(3,entity,step));
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity));
            set(traceHandle(entity),'YData', yData(:,entity));
            set(traceHandle(entity),'ZData', zData(:,entity));
        end
        % UPDATE TIMESTAMP ANNOTATION
        set(clockHandle,'String',sprintf('Time: %ss',num2str(SIM.TIME.timeVector(1,step))));
    end
    drawnow();
    
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    im = frame2im(F(frame));
    [imind,cm] = rgb2ind(im,256);
    
    % APPEND THE FRAMES TO GIF
    if frame == 1
        imwrite(imind,cm,filePath,'gif','Loopcount',inf,'DelayTime',(1/fps));
    else
        imwrite(imind,cm,filePath,'gif','WriteMode','append','DelayTime',(1/fps));
    end
    frame = frame + 1;                  % Move to next frame
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE 3D TRAJECTORY TRAIS AS A VIDEO [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricAvi(SIM,DATA,currentFigure)
% This function generates a video of the complete simulation trajectories
% over the complete series.

% OUTPUT FILE
fileName = strcat(SIM.outputPath,'isometricFigure','.avi');

% DATA CONTAINERS
% We need to build a matrix of states*tailLength*IDs,step
globalTraces = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),10,DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
globalMarker = NaN(10,SIM.totalObjects,SIM.TIME.numSteps); 

% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.endStep
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerState = objectStates(:,step);                            % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(:,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(:,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerState = objectStates(:,end);                             % Step position
            tailTrace = objectStates(:,((activeSteps + 1) - tailStepLength):end);
        end
                
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        globalMarker(:,indexValue,step) = markerState;      
        
        % BUILD MATICES OF TAIL COORDINATES
        % We need to build a matrix of dimensions
        % [tailLength*states*IDs,step], the resulting plot matrices must be
        % of dimensions [state(1)*tailLength,IDs]
        globalTraces(1:size(tailTrace,2),:,indexValue,step) = tailTrace';
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Timelapse');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');
set(figureHandle,'Visible','off');

% CONFIGURE VIDEO SETTINGS
fps = 50; % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.endTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end
% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(SIM.TIME.endStep/numFrames);                         % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end
vidFile = VideoWriter(fileName,'Motion JPEG AVI');
vidFile.FrameRate = fps;                                                   % Set the videos fps to match the sample rate
vidFile.Quality = 50;

open(vidFile);
frameSet = moviein(numFrames);                                             % Pre-allocate video object for frames
F(numFrames) = struct('cdata',[],'colormap',[]);

% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));
ax = gca(); 
ax.NextPlot = 'replaceChildren';

frame = 1; viewPoint = -45;
for step = 1:SIM.TIME.endStep
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % STEP DATA
    xData = squeeze(globalTraces(:,1,:,step));
    yData = squeeze(globalTraces(:,2,:,step));
    zData = squeeze(globalTraces(:,3,:,step));
    % BUILD THE FRAME 
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % INITIAL TRACE PLOT
        traceHandle = plot3(ax,xData,yData,zData);
        % UPDATE OBJECT HANDLES
        for ID1 = 1:SIM.totalObjects
            if isstruct(SIM.OBJECTS(ID1).patch)
               % GET THE GLOBAL POSE AT STEP
               [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,ID1,step));
               % BUILD MARKER
               markerHandle(ID1) = patch(ax,'Vertices',SIM.OBJECTS(ID1).patch.vertices*R2 + globalMarker(1:3,ID1,step)',...
                   'Faces',SIM.OBJECTS(ID1).patch.faces,...
                   'FaceColor',SIM.OBJECTS(ID1).colour,...
                   'EdgeColor',DATA.figureProperties.EdgeColor,...
                   'EdgeAlpha',DATA.figureProperties.EdgeAlpha,...  
                   'FaceLighting',DATA.figureProperties.FaceLighting,...
                   'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(ID1) = plot3(ax,globalMarker(1,ID1,step),globalMarker(2,ID1,step),globalMarker(3,ID1,step));
                % ASSIGN MARKER PROPERTIES
                markerHandle(ID1).Marker = SIM.OBJECTS(ID1).symbol;
                markerHandle(ID1).MarkerSize = DATA.figureProperties.MarkerSize;
                markerHandle(ID1).MarkerFaceColor = SIM.OBJECTS(ID1).colour;
                markerHandle(ID1).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
                markerHandle(ID1).Color = SIM.OBJECTS(ID1).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(ID1).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(ID1).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(ID1).Color = SIM.OBJECTS(ID1).colour;
            % ANNOTATION WITH THE CURRENT TIME
            clockHandle = annotation('textbox',[0.02 0.05 0.15 0.04],'String','Time:','FitBoxToText','off');
            set(clockHandle,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize,'FontSmoothing','on');
        xlabel('x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
        ylabel('y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
        zlabel('z(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
        set(ax,'Color',DATA.figureProperties.axesColor);
        set(ax,'GridAlpha',0.25,'GridColor','k');
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        set(ax,'outerposition',[0.02 0.1 1 0.88]);
        grid on; box on; hold off;
        legend(markerHandle,DATA.figureProperties.legendEntries,'Location','northeastoutside');
        view([-45 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            if isstruct(SIM.OBJECTS(entity).patch)
                % GET THE GLOBAL POSE AT STEP
                [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,entity,step));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + globalMarker(1:3,entity,step)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',globalMarker(1,entity,step));
                set(markerHandle(entity),'YData',globalMarker(2,entity,step));
                set(markerHandle(entity),'ZData',globalMarker(3,entity,step));
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity));
            set(traceHandle(entity),'YData', yData(:,entity));
            set(traceHandle(entity),'ZData', zData(:,entity));
        end
        % UPDATE TIMESTAMP ANNOTATION
        set(clockHandle,'String',sprintf('Time: %ss',num2str(SIM.TIME.timeVector(1,step))));
        view([viewPoint 50]);
    end
    drawnow();
    
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    writeVideo(vidFile,F(frame));       % Write frame to video
    % INCREMENT
    viewPoint = viewPoint + 0.05;
    frame = frame + 1;                  % Move to next frame
end
close(vidFile);
% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end

% //////////////////////////// TIMING FIGURES /////////////////////////////
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
figurePath = strcat(SIM.outputPath,'computationTimes');
figureHandle = figure('Name','OpenMAS Computation Time Series');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);

for ind = 1:SIM.totalObjects
    % GET OBJECT DATA
    entity = DATA.objectIndex{ind};
    % TEST FOR THE 'DATA' FIELD
    if ~isprop(entity,'DATA') || isempty(entity.DATA)
        continue
    end
    % GET THE EQUIVALENT SIM OBJECT
    IDvector = [SIM.OBJECTS.objectID];
    SIMobject = SIM.OBJECTS(IDvector == entity.objectID);
    if ~isfield(entity.DATA,'algorithm_dt')
        continue
    else
        % GET THE TIMESERIES DATA
        dt_timeSeries = entity.DATA.algorithm_dt*1E3; 
    end

    % LEGEND LABEL
    displayString = sprintf('[ID:%s] %s',num2str(SIMobject.objectID),SIMobject.name);
    % GENERATE FIGURES
    plot(DATA.timeVector(1:SIM.TIME.endStep),dt_timeSeries,...
         'Color',SIMobject.colour,...
         'LineWidth',DATA.figureProperties.LineWidth,...
         'DisplayName',displayString);
    hold on;
end
grid on; box on;
titleString = sprintf('Agent Computation times over a period of %ss',num2str(SIM.TIME.endTime));
title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize,'FontSmoothing','on');
xlabel('t (s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
ylabel('Computation Time(ms)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
legend('Location','northeastoutside');
set(gca,'FontSize',DATA.figureProperties.axisFontSize,'FontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
set(gca,'Color',DATA.figureProperties.axesColor);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;   
end

% ////////////////////////// MONTE-CARLO FIGURES //////////////////////////
% PLOT THE MEAN COMPUTATION TIME TIME-SERIES
function [figureHandle] = get_MonteCarloTimeSeries(cycleSample)
% This function generates the mean computation time - time
% series.

% BRING UP THE FIRST SIMPLE
[cycleSample] = obj.importCycleData(1);
exampleDATA = cycleSample.DATA;
exampleMETA = cycleSample.META;

% META FIGURE PROPERTIES
titleString = sprintf('Monte-Carlo (%s cycle) mean computation times',num2str(obj.cycles));
figurePath = strcat(obj.sessionPath,'\','computationTimeSeries');
figureHandle = figure('Name','Monte-Carlo Analysis: Agent mean computation times');
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.1, 0.90, 0.88]);
set(figureHandle,'Position',exampleDATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',exampleDATA.figureProperties.backGroundColor);          % Background colour

% GENERATE THE FIGURE
axesHandle = gca;
displayData = agentData.meanLoopTimeSeries*100;                 % Convert s to ms
lineHandle = plot(axesHandle,agentData.meanTimeVector,displayData);
set(lineHandle,'LineWidth',exampleDATA.figureProperties.LineWidth);

% FIGURE PROPERTIES
title(titleString,'fontweight',exampleDATA.figureProperties.fontWeight,'fontsize',exampleDATA.figureProperties.titleFontSize,'FontSmoothing','on');
xlabel(axesHandle,'t (s)','fontweight',exampleDATA.figureProperties.fontWeight,'fontSize',exampleDATA.figureProperties.axisFontSize,'FontSmoothing','on');
ylabel(axesHandle,'Computation Time (ms)','fontweight',exampleDATA.figureProperties.fontWeight,'fontSize',exampleDATA.figureProperties.axisFontSize,'FontSmoothing','on');
set(axesHandle,'FontSize',exampleDATA.figureProperties.axisFontSize,'fontWeight',exampleDATA.figureProperties.fontWeight,'FontSmoothing','on');
set(axesHandle,'GridAlpha',0.25,'GridColor','k');
xlim([ 0 exampleMETA.TIME.endTime]);
%            set(axesHandle,'outerposition',[0.01 0.05 1 0.88]);
grid on; box on; hold off;

% SAVE FIGURE TO OUTPUT DIRECTORY
savefig(figureHandle,figurePath);

% SAVE AS PDF
set(figureHandle,'Units','Inches');
pos = get(figureHandle,'Position');
set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
print(figureHandle,figurePath,'-dpdf','-r0');
end