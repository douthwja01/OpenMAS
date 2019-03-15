%% OPENMAS FIGURE GENERATOR (OMAS_figureGenerator.m) %%%%%%%%%%%%%%%%
% This script contains an index of the figures that can be requested from
% the simulation using the SIM.figures attribute.

% Author: James A. Douthwaite 23/11/2018

%% INDEX FUNCTION
function [figureNumber] = OMAS_figureGenerator(SIM,objectIndex,DATA,figureNumber,figureLabel)
% INPUTS:
% DATA         - The simulation output DATA structure
% figureNumber - The current figure number
% figureLabel  - The figure identifier
% OUTPUT:
% figureNumber - The updated figure number

% PARAMETERISE ALL THE FIGURES BY CURRENT CONFIGURATION
[DATA.figureProperties] = OMAS_figureProperties(SIM,DATA);
% GENERATE THE REDUCED STATE TRAJECTORY FILE
generateTrajectoryTempFile(SIM,DATA);
% DETERMINE WHICH FIGURE IS TO BE GENERATED
switch upper(char(figureLabel))
    case 'ALL'
        fprintf('[%s]\tAll output figures requested.\n',SIM.phase);
        % MOVE THROUGH THE COMPLETE FIGURE VECTOR
        figureVector = {'EVENTS','COLLISIONS','TRAJECTORIES','SEPARATIONS',...
                        'CLOSEST','INPUTS','PLAN','FIG','GIF','AVI','TIMES'};
        for fig = 1:length(figureVector)
            [figureNumber] = OMAS_figureGenerator(SIM,objectIndex,DATA,figureNumber,figureVector{fig});
        end
        close all;  % Kill figures
        fprintf('[%s]\n[%s]\tAll figures pushed to output directory: \n[%s]\t%s',...
                SIM.phase,SIM.phase,SIM.phase,SIM.outputPath);

    case 'EVENTS'
        fprintf('[%s]\tGenerating the event overview figure.\n',SIM.phase);
        [figureNumber,~] = get_eventOverview(SIM,DATA,figureNumber);  
        
    case 'COLLISIONS'
        fprintf('[%s]\tGenerating collision overview.\n',SIM.phase);
        figureSet = [];
        % CHECK COLLISIONS OCCURED
        if ~isfield(DATA,'uniqueCollisions') || isempty(DATA.uniqueCollisions)
            warning('[%s]\t...No collision data available.\n',SIM.phase);
            return
        end
        for collisionNumber = 1:numel(DATA.uniqueCollisions)
            [figureNumber,figureSet(collisionNumber)] = get_objectCollision(SIM,objectIndex,DATA,figureNumber,DATA.uniqueCollisions(collisionNumber));
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Collision Overview');
        set(windowHandle,'Position',DATA.figureProperties.windowSettings);        % Maximise the figure in the tab
        savefig(windowHandle,strcat(SIM.outputPath,'collision_overview'));        % Save the output figure
        
    case 'TRAJECTORIES'
        fprintf('[%s]\tGenerating global trajectory figure.\n',SIM.phase);
        figureSet = [];
        for objectNum = 1:DATA.totalObjects
            [figureNumber,figureSet(objectNum)] = get_objectTrajectory(SIM,DATA,figureNumber,objectNum);
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Trajectory Overview');
        set(windowHandle,'Position',DATA.figureProperties.windowSettings);        % Maximise the figure in the tab
        savefig(windowHandle,strcat(SIM.outputPath,'globalTrajectory_overview')); % Save the output figure
        
    case 'SEPARATIONS' 
        fprintf('[%s]\tGenerating object separation figure(s).\n',SIM.phase);
        figureSet = [];
        
        %if SIM.totalAgents < 1
        if (SIM.totalObjects - SIM.totalWaypoints) < 2
            warning('There must be at least two collidable objects to plot seperation data.\n');
            return
        end
        
        for objectNum = 1:(SIM.totalObjects - SIM.totalWaypoints)
            [figureNumber,figureSet(objectNum)] = get_objectSeparationFigure(SIM,DATA,figureNumber,objectNum);
        end
        
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Separation Overview');
        set(windowHandle,'Position', DATA.figureProperties.windowSettings);       % Maximise the figure in the tab
        savefig(windowHandle,strcat(SIM.outputPath,'separations_overview'));      % Save the output figure
        
    case 'CLOSEST'
        fprintf('[%s]\tGenerating closest separation figure.\n',SIM.phase);
        [figureNumber,~] = get_minimumSeparationFigure(SIM,DATA,figureNumber);
        
    case 'INPUTS'
        fprintf('[%s}\tGenerating control input figure.\n',SIM.phase);
        [figureNumber,~] = get_agentControlInputs(SIM,objectIndex,DATA,figureNumber);
        
    case 'PLAN'
        fprintf('[%s}\tGenerating top-down 2D(plan) figure.\n',SIM.phase);
        [figureNumber,~] = get_topDownView(SIM,objectIndex,DATA,figureNumber);
        
    case 'FIG' 
        fprintf('[%s]\tGenerating isometric trajectory figure.\n',SIM.phase);
        [figureNumber,~] = get_isometricFigure(SIM,objectIndex,DATA,figureNumber);

    case 'GIF'
        fprintf('[%s]\tGenerating trajectory gif file.\n',SIM.phase);
        [figureNumber,~] = get_isometricGif(SIM,objectIndex,DATA,figureNumber);
        
    case 'AVI'
        fprintf('[%s]\tGenerating trajectory avi file.\n',SIM.phase);
        [figureNumber,~] = get_isometricAvi(SIM,objectIndex,DATA,figureNumber);
        
    case 'TIMES'
        fprintf('[%s]\tGenerating computation timeseries figure.\n',SIM.phase);
        [figureNumber,~] = get_computationTimes(SIM,objectIndex,DATA,figureNumber);
        
    case 'NONE'
        fprintf('[%s]\tFigure output supressed.\n',SIM.phase);
        
    otherwise
        warning('[WARNING] Ignoring figure request "%s": Figure unknown.',char(figureLabel));
end
end

% PREPARATION FUNCTIONS 
function [filePath] = generateTrajectoryTempFile(SIM,DATA)
% This function generates/appends to a temp file containing the trajectory 
% states to be plotted in annimations. The number of states are reduced
% according to the desired frame rate and are lower resolution than the
% true number of steps in some cases.
% OUTPUT:
% filePath        - The path to the system file
% objectIDx       - The trajectory data reduced to the desired number of frames.
% frameTimeVector - The reduced time vector associated with the frames.

% CREATE A TEMP FILE FOR THE GENERATION OF TRAJECTORY ELEMENTS
numRep = floor(SIM.TIME.endStep/DATA.figureProperties.stepsPerFrame);          % Handle numsteps not divisable by stepsPerFrame
logicalMatrix = logical(zeros(1,numRep*DATA.figureProperties.stepsPerFrame));  % Frame selection matrix

% FOR EACH OBJECT, GET THE STATE SET
outputStructure = struct();
for objectNo = 1:DATA.totalObjects
    % GET THE COMPLETE STATE SET
    completeStateSet = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,...
                                                  SIM.globalIDvector,...
                                                  SIM.OBJECTS(objectNo).objectID,...
                                                  inf);                    % All valid states for the object
    completeStateSet = completeStateSet(:,1:numRep*DATA.figureProperties.stepsPerFrame);
    
    % GET THE INDICES OF EACH FRAME IN THE STATE SPACE
    for num = 1:numRep
        logicalMatrix(1,num*DATA.figureProperties.stepsPerFrame) = logical(true);
    end
    % SELECT THE STATES FROM THE COMPLETE SET
    statesToPlot = completeStateSet(:,logicalMatrix);
    % GENERATE DATA STRUCTURE TO BE SAVED
    outputStructure.(sprintf('objectID%d',SIM.OBJECTS(objectNo).objectID)) = statesToPlot;
end
% GET THE TIME VECTOR FOR THE FRAME SET
frameTimeVector = SIM.TIME.timeVector(1:numRep*DATA.figureProperties.stepsPerFrame);
outputStructure.frameTimeVector = frameTimeVector(logicalMatrix);

% SAVE A TEMPORARY FILE WITH THE STATES OF EACH OBJECT
filePath = [SIM.outputPath,SIM.systemFile];
save(filePath,'-struct','outputStructure','-append');
end

%% /////////////////// COLLISION DATA FIGURE GENERATION ///////////////////
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
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-Label
xlabel(h,...
    't (s)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
xlim([-SIM.TIME.dt,SIM.TIME.endTime])
% Y-Label
ylabel(h,...
    'Number of Events',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Axes
set(h,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor);
% Legend
legend(h,eventFields,...
    'Location','northeastoutside',...
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
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end
% PLOT THE COLLISION EVENT DESCRIPTION FIGURES
function [currentFigure,figureHandle] = get_objectCollision(SIM,objectIndex,DATA,currentFigure,collisionEvent)
% This function is designed to move through the collsion event history and
% generate figure for each of the collisions.

% Generate the figure
figureHandle = figure('Name',sprintf('[t=%.0fs] %s & %s',collisionEvent.time,collisionEvent.name_A,collisionEvent.name_B));
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour
ax = axes(figureHandle);
hold on;
view([70 25]);
titlestr = sprintf('Collision Event between %s and %s at t=%s',...         % Declare title string for figure
           collisionEvent.name_A,collisionEvent.name_B,num2str(collisionEvent.time));

% /////////////////////////// ENTITY DATA /////////////////////////////////
% EVENT, OBJECT AND META DATA
META_A = SIM.OBJECTS(SIM.globalIDvector == collisionEvent.objectID_A);
META_B = SIM.OBJECTS(SIM.globalIDvector == collisionEvent.objectID_B);     % The associated META datas
object_A = objectIndex{SIM.globalIDvector == collisionEvent.objectID_A};
object_B = objectIndex{SIM.globalIDvector == collisionEvent.objectID_B};   % The associated object structures
geometry_A = object_A.GEOMETRY;
geometry_B = object_B.GEOMETRY;
[objectAStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,META_A.objectID,inf);
[objectBStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,META_B.objectID,inf);

% We need to plot:
% - The object geometries
% - The object hitboxes

% ////////////////////// PLOT THE OBJECT GEOMETRY /////////////////////////
% SUBSTITUTE A'S GEOMETRY IF NECESSARY
if size(object_A.GEOMETRY.vertices,1) < 1
    [geometry_A] = OMAS_graphics.defineSphere(zeros(3,1),META_A.radius,10);
end
% GET A'S ROTATIONS
R_A = OMAS_geometry.quaternionToRotationMatrix(collisionEvent.state_A(7:10)); 
% OBJECT A'S GEOMETRY
patch(ax,...
    'Vertices',(geometry_A.vertices*R_A + collisionEvent.state_A(1:3)'),...
    'Faces',geometry_A.faces,...
    'FaceColor',META_A.colour,...
    'EdgeColor','k',...
    'EdgeAlpha',DATA.figureProperties.edgeAlpha,...
    'FaceAlpha',DATA.figureProperties.faceAlpha,...
    'FaceLighting','gouraud',...
    'LineWidth',DATA.figureProperties.lineWidth);

% SUBSTITUTE B'S GEOMETRY IF NECESSARY
if size(object_B.GEOMETRY.vertices,1) < 1
    [geometry_B] = OMAS_graphics.defineSphere(zeros(3,1),META_B.radius,10); 
end
% GET B'S ROTATIONS
R_B = OMAS_geometry.quaternionToRotationMatrix(collisionEvent.state_B(7:10)); 
% OBJECT B'S GEOMETRY
patch(ax,...
    'Vertices',(geometry_B.vertices*R_B + collisionEvent.state_B(1:3)'),...
    'Faces',geometry_B.faces,...
    'FaceColor',META_B.colour,...
    'EdgeColor','k',...
    'EdgeAlpha',DATA.figureProperties.edgeAlpha,...
    'FaceAlpha',DATA.figureProperties.faceAlpha,...
    'FaceLighting','gouraud',...
    'LineWidth',DATA.figureProperties.lineWidth);

% ////////////////////// PLOT THE OBJECT'S HITBOX /////////////////////////
% Get the hit-box geometry
[hitBoxGeometry_A] = OMAS_graphics.getHitBoxGeometry(object_A.VIRTUAL,object_A.GEOMETRY);
% DEFINE A's HITBOX
patch(ax,...
    'Vertices',(hitBoxGeometry_A.vertices + collisionEvent.state_A(1:3)'),...
    'Faces',hitBoxGeometry_A.faces,...
    'FaceColor','r',...
    'EdgeColor','k',...
    'EdgeAlpha',DATA.figureProperties.edgeAlpha*0.8,...
    'FaceAlpha',DATA.figureProperties.faceAlpha*0.8,...
    'FaceLighting','gouraud',...
    'LineWidth',DATA.figureProperties.lineWidth);

% Get the hit-box geometry
[hitBoxGeometry_B] = OMAS_graphics.getHitBoxGeometry(object_B.VIRTUAL,object_B.GEOMETRY);
% DEFINE B's HITBOX
patch(ax,...
    'Vertices',(hitBoxGeometry_B.vertices + collisionEvent.state_B(1:3)'),...
    'Faces',hitBoxGeometry_B.faces,...
    'FaceColor','r',...
    'EdgeColor','k',...
    'EdgeAlpha',DATA.figureProperties.edgeAlpha*0.8,...
    'FaceAlpha',DATA.figureProperties.faceAlpha*0.8,...
    'FaceLighting','gouraud',...
    'LineWidth',DATA.figureProperties.lineWidth);

% ////////////////////// PLOT THE TRAJECTORY TRAILS ///////////////////////
plot3(ax,objectAStates(1,:),objectAStates(2,:),objectAStates(3,:),...
    'LineStyle',DATA.figureProperties.lineStyle,...
    'LineWidth',DATA.figureProperties.lineWidth,...
    'Color',META_A.colour);
plot3(ax,objectBStates(1,:),objectBStates(2,:),objectBStates(3,:),...
    'LineStyle',DATA.figureProperties.lineStyle,...
    'LineWidth',DATA.figureProperties.lineWidth,...
    'Color',META_B.colour);

% ////////////////////// PLOT THE VELOCITY VECTORS ////////////////////////
% PLOT THE OBJECT VELOCITY VECTORS
q = quiver3(collisionEvent.state_A(1),...
    collisionEvent.state_A(2),...
    collisionEvent.state_A(3),...
    collisionEvent.state_A(4),...
    collisionEvent.state_A(5),...
    collisionEvent.state_A(6),'r');
q.AutoScaleFactor = 1;
q.LineWidth = DATA.figureProperties.lineWidth;
q = quiver3(collisionEvent.state_B(1),...
    collisionEvent.state_B(2),...
    collisionEvent.state_B(3),...
    collisionEvent.state_B(4),...
    collisionEvent.state_B(5),...
    collisionEvent.state_B(6),'g');
q.AutoScaleFactor = 1;
q.LineWidth = DATA.figureProperties.lineWidth;

% ///////////////////////// ADD LABEL ANNOTATIONS /////////////////////////
annotationText = sprintf('    %s [ID-%d]',META_A.name,META_A.objectID);
text(collisionEvent.state_A(1),collisionEvent.state_A(2),collisionEvent.state_A(3),annotationText);
annotationText = sprintf('    %s [ID-%d]',META_B.name,META_B.objectID);
text(collisionEvent.state_B(1),collisionEvent.state_B(2),collisionEvent.state_B(3),annotationText);

% DEFINE THE PLOT LAYOUT
% Title
title(ax,titlestr,...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-label
xlabel(ax,'x(m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Y-Label
ylabel(ax,'y(m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Z-Label
zlabel(ax,'z(m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Axes
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
axisLimit = (DATA.figureProperties.objectMaximalRadii + DATA.figureProperties.maxAbsPosition);
xlim(ax,[-axisLimit,axisLimit]);
ylim(ax,[-axisLimit,axisLimit]);
zlim(ax,[-axisLimit,axisLimit]);
axis vis3d equal;    
grid on;  box on;
hold off;

% SAVE THE OUTPUT FIGURE
filename = sprintf('collisionEvent_[ID-%d] %s_[ID-%d] %s.fig',...
                    collisionEvent.objectID_A,collisionEvent.name_A,...
                    collisionEvent.objectID_B,collisionEvent.name_B);
savefig(figureHandle,strcat(SIM.outputPath,filename));
% ITERATE PLOT
currentFigure = currentFigure + 1;
end

%% ////////////////// TRAJECTORY DATA FIGURE GENERATION ///////////////////
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
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1,0.1,0.85,0.82]);     

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
[objectStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(objectNum).objectID,inf);
% STATE NAME VECTOR
stateTags = {{'$x$','$(m)$'},{'$y$','$(m)$'},{'$z$','$(m)$'},...
             {'$u$','$(m/s)$'},{'$v$','$(m/s)$'},{'$w$','$(m/s)$'},...
             {'$\phi$','$(rad)$'},{'$\theta$','$(rad)$'},{'$\psi$','$(rad)$'}};
% CONVERT QUATERNION STATE TO EULER ANGLE STATES         
eulerStates = zeros(9,size(objectStates,2));         
for i = 1:size(objectStates,2)  
    eulerStates(1:6,i) = objectStates(1:6,i);
    eulerStates(7:9,i) = OMAS_geometry.quaternionToEulers(objectStates(7:10,i));         
end            

clear objectStates;

% FOR EACH OBJECTS PERSPECTIVE
stateVectorLength = size(eulerStates,1);
for n = 1:stateVectorLength       
    % CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
    plotCellB = n*plotCellWidth;                                    % The end of the plot
    plotLocation = subplot(stateVectorLength,plotCellWidth,[plotCellA plotCellB]);
    
    % PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,DATA.timeVector(1:SIM.TIME.endStep),...
         eulerStates(n,1:SIM.TIME.endStep),...
         'LineStyle','-',...
         'LineWidth',DATA.figureProperties.lineWidth,...
         'Color','b');
    % Y-axes
    ylabel(plotLocation,stateTags{n},...
        'Interpreter',DATA.figureProperties.interpreter,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on');
    % Axes
    set(plotLocation,...
        'TickLabelInterpreter',DATA.figureProperties.interpreter,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on',...
        'Color',DATA.figureProperties.axesColor);
    grid on; box on; grid minor;
    plotLocation.XAxis.MinorTickValues = plotLocation.XAxis.Limits(1):SIM.TIME.dt:plotLocation.XAxis.Limits(2);
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if n == 1
        % Title
        title(plotLocation,titlestr,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'Fontsize',DATA.figureProperties.titleFontSize,...
            'FontSmoothing','on');   
        % Append title to first subplot
    end
    % Prevent overlap of x-label
    if n == stateVectorLength
        % X-axis
        xlabel(plotLocation,'t (s)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
    else
        set(plotLocation,'XTickLabel',[]); 
    end
    % Move to next subplot location 
    plotCellA = plotCellA + plotCellWidth;
end
hold off; 

% SAVE THE OUTPUT FIGURE
fileName = strcat(SIM.outputPath,sprintf('globalTrajectory_[ID-%.0f] %s',SIM.OBJECTS(objectNum).objectID,SIM.OBJECTS(objectNum).name));
set(figureHandle,'Visible','on');                                        % Make it visable for saving
savefig(figureHandle,fileName);                                            % As matlab figure 
set(figureHandle,'Visible','off');
% ITERATE PLOT
currentFigure = currentFigure + 1;
end
% GET THE AGENT SEPERATION SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = get_objectSeparationFigure(SIM,DATA,currentFigure,objectNum)
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
% GET THE CLOSEST PROXIMITY FIGURES
function [currentFigure,figureHandle] = get_minimumSeparationFigure(SIM,DATA,currentFigure)
% This function generates a figure with the smallest separations observed
% in the simulation and plots their separation. This plot only plots the
% separation between objects that have a form of hitbo (i.e. ~none).

% INPUT HANDLING
figureHandle = [];
if ~any([SIM.OBJECTS(:).hitBox] ~= OMAS_hitBoxType.none) 
    warning('There must be at least two collidable objects to plot seperation data.');
    return
end
if ~any([SIM.OBJECTS(:).type] ~= OMAS_objectType.waypoint) 
    warning('There are no objects that are not waypoints.');
    return
end

% Get the META structures for all hitbox objects
hitBoxOBJECTS = SIM.OBJECTS([SIM.OBJECTS.hitBox] ~= OMAS_hitBoxType.none);   % The objects with valid hitboxes
hitBoxOBJECTS = hitBoxOBJECTS([hitBoxOBJECTS.type] ~= OMAS_objectType.waypoint);
% Get objects that are valid for consideration
hitBoxIDvector = [hitBoxOBJECTS.objectID];
% Output container
separationTimeSeries = inf(numel(hitBoxIDvector),numel(hitBoxIDvector),SIM.TIME.endStep);
for i = 1:numel(hitBoxOBJECTS)
    % Get the state trace for object A
    objectStatesA = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,uint8(hitBoxIDvector(i)),inf);
    objectStatesA = objectStatesA(:,1:SIM.TIME.endStep);
    for j = 1:numel(hitBoxOBJECTS)
        if hitBoxIDvector(i) == hitBoxIDvector(j)
            continue
        end
        % Get the state trace for object B
        objectStatesB = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,uint8(hitBoxIDvector(j)),inf);
        objectStatesB = objectStatesB(:,1:SIM.TIME.endStep);
        % Calculate the separation across the time period
        centroidSeparations = objectStatesB(1:3,:) - objectStatesA(1:3,:);
        separationTimeSeries(i,j,:) = sqrt(sum(centroidSeparations.^2,1));
    end
end

% GET THE MINIMUM SEPERATIONS FOR EACH AGENT
% If the number of obstacle/agents exceeds 10 objects, we are only
% interested in the 10 interactions that came the closest to collision.
maxDisplayObjects = 10;

% MIN AND MAXIMUM SEPERATION MATRICES [IDA by IDB]
minABMatrix = min(separationTimeSeries,[],3);       % Minimum separations over the timeseries
maxABMatrix = max(separationTimeSeries,[],3);       % Maximum separations over the timeseries
maxABMatrix(isinf(maxABMatrix)) = NaN;              % Remove the 'inf' cases
[h_maxVals,~] = max(maxABMatrix,[],2);              % The minimum seperations, horezontal indicies
maxProximity  = max(h_maxVals);                     % Maximal seperation
% Find the indices 
[~,indexOfClosestObject] = min(minABMatrix,[],2);   % The minimum with respect to the IDvector

% Define a vector of comparative IDs
interactionsByID = [hitBoxIDvector',hitBoxIDvector(indexOfClosestObject)'];
interactionsByID = sort(interactionsByID,2);
uniqueIDInteractions = unique(interactionsByID,'rows');

% Limit display objects
if size(uniqueIDInteractions,1) > maxDisplayObjects
    uniqueIDInteractions = uniqueIDInteractions(1:maxDisplayObjects,:);
end

% DEFINE THE FIGURE
figurePath = strcat(SIM.outputPath,'minimumSeparations');   
figureHandle = figure('Name','OpenMAS point of closest approach'); 
setappdata(figureHandle,'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
ax = axes(figureHandle);                                                   % Begin generating the figure
% ADD FINAL PLOT ATTRIBUTES
% Title
titleString = sprintf('Global minimum separations over a period of %ss',num2str(SIM.TIME.endTime));
title(ax,titleString,...
    'Interpreter',DATA.figureProperties.interpreter,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSmoothing','On');
% X-axes
xlabel(ax,'t (s)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSmoothing','On');
% Y-axes
ylabel(ax,'Separation (m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSmoothing','On');
% Axes properties
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight,...
    'Color',DATA.figureProperties.axesColor,...
    'FontSmoothing','On',...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
ylim(ax,[0 1.1*maxProximity]);                                             % Use the maximum recorded separation 
grid on; box on; hold on; grid minor;
% Show the timestep difference in the figure
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):SIM.TIME.dt:ax.XAxis.Limits(2);

% Each of the uniqe interactions
legendEntries = {}; legendCounter = 1;
for i = 1:size(uniqueIDInteractions,1)
    % The associate META objects
    META_A = hitBoxOBJECTS([hitBoxOBJECTS.objectID] == uniqueIDInteractions(i,1));
    META_B = hitBoxOBJECTS([hitBoxOBJECTS.objectID] == uniqueIDInteractions(i,2));
    ind_A  = find(hitBoxIDvector == META_A.objectID);
    ind_B  = find(hitBoxIDvector == META_B.objectID);
    % Get the separation tracing 
    separationTrace = squeeze(separationTimeSeries(ind_A,ind_B,:)); 
    % Collective radii
    criticalLimit = META_A.radius + META_B.radius;
    % Build the labels
    label1 = sprintf('%s[ID-%d]',META_A.name,META_A.objectID);
    label2 = sprintf('%s[ID-%d]',META_B.name,META_B.objectID);
    % Create the legend entries
    legendEntries{legendCounter} = strcat(label1,' - ',label2);
    % PLOT THE SEPERATIONS
    plot(ax,DATA.timeVector(:,1:SIM.TIME.endStep),...
         separationTrace',...
        'Color',META_A.colour,...
        'LineStyle','-',...
        'LineWidth',DATA.figureProperties.lineWidth);
    % INCREMENT THE LEGEND COUNTER
    legendCounter = legendCounter + 1;                             % Increment legend entry
end

% THE COLLISION CONDITION
refHandle = refline(ax,0,criticalLimit);                                % Adds a reference line with slope m and intercept b to the current axes.
set(refHandle,'color','k','LineStyle','--','LineWidth',DATA.figureProperties.lineWidth);

% Add legend entries
legendEntries{legendCounter} = 'Collision Boundary';
legend(ax,...
    legendEntries,...
    'Location','northEast',...
    'Interpreter',DATA.figureProperties.interpreter);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;    
end
% AGENT CONTROL INPUT TRAJECTORIES
function [currentFigure,figureHandle] = get_agentControlInputs(SIM,objectIndex,DATA,currentFigure)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

figureHandle = [];

% SANITY CHECK 1
if DATA.totalAgents == 0
    warning('No agent data available.'); 
    return 
else
    % LOGICALLY SELECT THE AGENTS FROM THE OBJECT INDEX 
    agentObjectIndex = objectIndex([SIM.OBJECTS.type] == OMAS_objectType.agent); % The agents themselves  
end

iter = 0;
for ID1 = 1:SIM.totalAgents
    % TRY TO GET THE REQUIRED PROPERTIES
    try
       testA = agentObjectIndex{ID1}.DATA.inputNames;
       testB = agentObjectIndex{ID1}.DATA.inputs;
    catch
       iter = iter + 1; 
    end
end

% SANITY CHECK 2
if iter == SIM.totalAgents
    warning('No agents with defined input timeseries. Use the "DATA.inputNames" and "DATA.inputs" fields if desired.');
    return 
end

% Declare title string for figure  
titlestr = sprintf('Agent control trajectories over a period of %ss',num2str(SIM.TIME.endTime));
figurePath = strcat(SIM.outputPath,'inputs');

% FIGURE META PROPERTIES
figureHandle = figure('Name','OpenMAS control inputs');
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.88, 0.85]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
plotCellA = 1; plotCellWidth = 4;                                          % The width of each figure, the start of the plot

% FOR EACH AGENT'S PERSPECTIVE
for ID1 = 1:SIM.totalAgents
    % GET EVALUATION OBJECT
    evalAgent = agentObjectIndex{ID1};
    % CHECK SKIP CONDITION
    if isempty(evalAgent.DATA)
        warning('[OUTPUT]\tProperty .DATA is empty for agent %s',evalAgent.name);
        continue
    end
    
    % GET TH AGENT-LOCAL 'DATA' STRUCTURE
    objectDATA = evalAgent.DATA;
    inputTrajectories = objectDATA.inputs;
    inputCount = size(inputTrajectories,1);
    
    % BUILD THE SUB-PLOT
    plotCellB = double(ID1*plotCellWidth);                                         % The end of the plot
    plotLocation = subplot(double(DATA.totalAgents),plotCellWidth,[plotCellA plotCellB]);
    set(plotLocation,'Color',DATA.figureProperties.axesColor);
    set(plotLocation,...
        'GridLineStyle','--',...
        'GridAlpha',0.25,...
        'GridColor','k');
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
    
    % PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,...
         DATA.timeVector(:,1:SIM.TIME.endStep),...
         inputTrajectories(:,1:SIM.TIME.endStep),...
        'LineStyle','-',...
        'LineWidth',DATA.figureProperties.lineWidth);
    % Y-Label
    ylabel(plotLocation,...
        sprintf('Inputs \n [ID-%s]',num2str(evalAgent.objectID)),...
        'Interpreter',DATA.figureProperties.interpreter,...
        'fontweight',DATA.figureProperties.fontWeight,...
        'fontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','On');
    % Axes
    set(plotLocation,...
        'TickLabelInterpreter',DATA.figureProperties.interpreter,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on',...
        'Color',DATA.figureProperties.axesColor,...
        'GridLineStyle','--',...
        'GridAlpha',0.25,...
        'GridColor','k');
    % Legend
    legend(legendEntries,...
        'Location','eastoutside',...
        'Interpreter',DATA.figureProperties.interpreter);
    
    grid on; box on; grid minor;
    % Show the timestep difference in the figure
    plotLocation.XAxis.MinorTickValues = plotLocation.XAxis.Limits(1):SIM.TIME.dt:plotLocation.XAxis.Limits(2);
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        % Title
        title(plotLocation,titlestr,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontsize',DATA.figureProperties.titleFontSize);                                          % Append title to first subplot
    end
    % Prevent overlap of x-label
    if ID1 == DATA.totalAgents
        % X-Label
        xlabel(plotLocation,'t (s)',...
            'Interpreter', DATA.figureProperties.interpreter,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','On');                                        % Append title to first subplot
    else
        set(plotLocation,'XTickLabel',[]); 
    end

    % Move to next subplot location       
    plotCellA = plotCellA + plotCellWidth; 
end
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end 
% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;  
end

%% ////////////////////// GENERAL OVERVIEW FIGURES ////////////////////////
% GET THE STANDARD 2D TRAJECTORY DATA.figureProperties [ UPDATED ]
function [currentFigure,figureHandle] = get_topDownView(SIM,objectIndex,DATA,currentFigure)

% FIGURE TITLE
figurePath = strcat(SIM.outputPath,'2Dtrajectories');

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS top-down image');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
ax = axes(figureHandle);
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW

legendCounter = 1; legendEntries = cell(DATA.totalObjects,1);

hold on; grid on; box on; grid minor;
for ID1 = 1:DATA.totalObjects    
    % GET OBJECT OVERVIEW DATA
    legendString = sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(ID1).objectID),SIM.OBJECTS(ID1).name);
    legendEntries(legendCounter) = {legendString};
    % THE OBJECT HANDLE
    objectHandle = objectIndex{SIM.OBJECTS(ID1).objectID == SIM.globalIDvector};
    % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
    [finalStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,SIM.TIME.endStep);
    finalPosition = finalStates(1:3,:);
    % LOCAL FIXED TO GLOBAL ROTATED
    R_final = OMAS_geometry.quaternionToRotationMatrix(finalStates(7:10)); 
    % DISPLAY THE OBJECT
    if numel(objectHandle.GEOMETRY.vertices) > 0
        patch(ax,'Vertices',objectHandle.GEOMETRY.vertices*R_final + finalPosition',...
            'Faces',objectHandle.GEOMETRY.faces,...
            'FaceColor',SIM.OBJECTS(ID1).colour,...
            'EdgeColor',DATA.figureProperties.edgeColor,...
            'EdgeAlpha',DATA.figureProperties.edgeAlpha,...  
            'FaceLighting',DATA.figureProperties.faceLighting,...
            'FaceAlpha',DATA.figureProperties.faceAlpha,...
            'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties
    else
        % PLOT THE TERMINAL POSITIONS
        plot3(ax,finalPosition(1),finalPosition(2),finalPosition(3),...
              'Marker',SIM.OBJECTS(ID1).symbol,...
              'MarkerSize',DATA.figureProperties.markerSize,...
              'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
              'MarkerEdgeColor',DATA.figureProperties.markerEdgeColor,...
              'LineWidth',DATA.figureProperties.lineWidth,...
              'LineStyle',DATA.figureProperties.lineStyle,...
              'Color',SIM.OBJECTS(ID1).colour); 
    end  
    legendCounter = legendCounter + 1;
end
legend('off')
hold on;
for ID1 = 1:DATA.totalObjects 
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    idleFlag = NaN('double');
    [objectStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,idleFlag);
    positions = objectStates(1:3,:);
    plot3(positions(1,:),positions(2,:),positions(3,:),...
          'LineStyle',DATA.figureProperties.lineStyle,...
          'LineWidth',DATA.figureProperties.lineWidth,...
          'Color',SIM.OBJECTS(ID1).colour);
      
    % ADD THE INTIAL POSITION REFERENCE POSITION 
    if SIM.OBJECTS(ID1).type == OMAS_objectType.agent
        initialPlotPosition = positions(1:2,1) - SIM.OBJECTS(ID1).radius; % The first position of the object
        rectangle('Position',[initialPlotPosition(1) initialPlotPosition(2) (2*SIM.OBJECTS(ID1).radius) (2*SIM.OBJECTS(ID1).radius)],...
                  'Curvature',[1 1],...
                  'FaceColor',SIM.OBJECTS(ID1).colour,...
                  'EdgeColor',DATA.figureProperties.edgeColor,...
                  'LineWidth',DATA.figureProperties.patchLineWidth)  
    end  
end

% Title
title(ax,...
    sprintf('Object trajectories over a period of %ss',num2str(SIM.TIME.endTime)),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-Label
xlabel(ax,'x(m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Y-Label
ylabel(ax,'y(m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Axes
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
% Legend
legend(legendEntries,...
    'location','northeastoutside',...
    'interpreter',DATA.figureProperties.interpreter);

% axisMins = DATA.figureProperties.axisMinimums(1:3) - DATA.figureProperties.objectMaximalRadii;
% axisMaxs = DATA.figureProperties.axisMaximums(1:3) + DATA.figureProperties.objectMaximalRadii;
% xlim(ax,[axisMins(1),axisMaxs(1)]);
% ylim(ax,[axisMins(2),axisMaxs(2)]);
% zlim(ax,[axisMins(3),axisMaxs(3)]);
axis square equal;

% axis square;     
view([0 90]);
% set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
% GET THE STANDARD 3D TRAJECTORY DATA.figureProperties [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricFigure(SIM,objectIndex,DATA,currentFigure)

% FIGURE TITLE
figurePath = strcat(SIM.outputPath,'isometricFigure');

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS isometric view');
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]); % MAXIMISE GRAPH SIZE IN WINDOW
ax = axes(figureHandle);
legendCounter = 1; legendEntries = cell(DATA.totalObjects,1);

hold on; grid on; box on; grid minor;
for ID1 = 1:DATA.totalObjects    
    % GET OBJECT OVERVIEW DATA
    legendString = sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(ID1).objectID),SIM.OBJECTS(ID1).name);
    legendEntries(legendCounter) = {legendString};
    % THE ASSOCIATED LOGIC
    objectID1 = objectIndex{SIM.globalIDvector == SIM.OBJECTS(ID1).objectID};
    % EXTRACT FINAL POSITION DATA FROM THE TRAJECTORY MATRIX
    [finalStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,SIM.TIME.endStep);
    finalPosition = finalStates(1:3,:);
    
    % LOCAL FIXED TO GLOBAL ROTATED
    R_final = OMAS_geometry.quaternionToRotationMatrix(finalStates(7:10));
    
    % DISPLAY THE OBJECT
    if numel(objectID1.GEOMETRY.vertices) > 0
        patch(ax,'Vertices',objectID1.GEOMETRY.vertices*R_final + finalPosition',...
            'Faces',objectID1.GEOMETRY.faces,...
            'FaceColor',SIM.OBJECTS(ID1).colour,...
            'EdgeColor',DATA.figureProperties.edgeColor,...
            'EdgeAlpha',DATA.figureProperties.edgeAlpha,...  
            'FaceLighting',DATA.figureProperties.faceLighting,...
            'FaceAlpha',DATA.figureProperties.faceAlpha,...
            'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties            % Patch properties
    else
        % PLOT THE TERMINAL POSITIONS
        plot3(ax,finalPosition(1),finalPosition(2),finalPosition(3),...
              'Marker',SIM.OBJECTS(ID1).symbol,...
              'MarkerSize',DATA.figureProperties.markerSize,...
              'MarkerFaceColor',SIM.OBJECTS(ID1).colour,...
              'MarkerEdgeColor',DATA.figureProperties.markerEdgeColor,...
              'LineWidth',DATA.figureProperties.lineWidth,...
              'LineStyle',DATA.figureProperties.lineStyle,...
              'Color',SIM.OBJECTS(ID1).colour); 
    end  
    legendCounter = legendCounter + 1;
end
legend('off')
hold on;
for ID1 = 1:DATA.totalObjects 
    % EXTRACT STATE TIME-SERIES DATA UPTO THE IDLE POINT
    idleFlag = NaN('double');
    [objectStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(ID1).objectID,idleFlag);
    positions = objectStates(1:3,:);
    plot3(positions(1,:),positions(2,:),positions(3,:),...
          'LineStyle',DATA.figureProperties.lineStyle,...
          'LineWidth',DATA.figureProperties.lineWidth,...
          'Color',SIM.OBJECTS(ID1).colour);
end

% Title
title(ax,...
    sprintf('Object trajectories over a period of %ss',num2str(SIM.TIME.endTime)),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-Label
xlabel(ax,'x (m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Y-Label
ylabel(ax,'y (m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Z-Label
zlabel(ax,'z (m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Axes
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
% Legend
legend(legendEntries,...
    'location','northeastoutside',...
    'interpreter',DATA.figureProperties.interpreter);

% xlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
% ylim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%zlim(ax,[-DATA.figureProperties.maxAbsPosition,DATA.figureProperties.maxAbsPosition]);
%xlim(ax,[DATA.figureProperties.axisMinimums(1)-0.1,DATA.figureProperties.axisMaximums(1)+0.1]);
%ylim(ax,[DATA.figureProperties.axisMinimums(2)-0.1,DATA.figureProperties.axisMaximums(2)+0.1]);
% zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);

axis vis3d equal;     
view([-24 36]);
set(ax,'outerposition',[0.05 0.15 1 0.68]);                               % Set the axes offset position in the figure window
grid on; grid minor;
hold off;

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

% FIGURE COMPLETE
currentFigure = currentFigure + 1;
% CLEAN UP
clearvars -except currentFigure figureHandle
end
% GET THE 3D TRAJECTORY TRAILS AS A GIF  [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricGif(SIM,objectIndex,DATA,currentFigure)
% This function generates an animated .gif representing the object
% trajectories over the complete timeseries.

% CONFIGURE THE PLOT ATTRIBUTES
filePath = strcat(SIM.outputPath,'isometricFigure.gif');
figureHandle = figure('Name','OpenMAS isometric timelapse (GIF)');
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));

% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position',DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
if DATA.figureProperties.publish
    set(figureHandle,'MenuBar','none');
    set(figureHandle,'ToolBar','none');
    set(figureHandle,'Visible','off');
end
% AXES SETTINGS
ax = axes(figureHandle);

% GET THE SIZE OF THE PLOTTED SET
% The data has already been pre-formatted so that the number of states
% aligns with the number of expected frames in the annimation.
% IMPORT THE OBJECT DATA FROM TEMP FILE
load([SIM.outputPath,SIM.systemFile]);
framesToPlot = size(eval(sprintf('objectID%d',SIM.OBJECTS(1).objectID)),2);   % The variable name in the workspace
F(framesToPlot) = struct('cdata',[],'colormap',[]);

% ax.NextPlot = 'replaceChildren';
inclinationAngle = 62;
viewPoint = -10;
viewRate = 0.05;
for frame = 1:framesToPlot
    hold on;
    % Must be done on a step by step base to get the annimations correct.
    
    % BUILD THE FRAME 
    if frame == 1
        % UPDATE OBJECT HANDLES
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT HANDLE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            globalStates = eval(variableLabel);                              % Evalutate object states
            
            % TRACE INITIAL POSITIONS
            traceHandle(entity) = plot3(ax,globalStates(1,frame),globalStates(2,frame),globalStates(3,frame));
            if numel(objectHandle.GEOMETRY.vertices) > 0
                % GET THE GLOBAL POSE AT STEP
                [R_frame] = OMAS_geometry.quaternionToRotationMatrix(globalStates(7:end,frame));
                % BUILD MARKER
                markerHandle(entity) = patch(ax,...
                    'Vertices',objectHandle.GEOMETRY.vertices*R_frame + globalStates(1:3,frame)',...
                    'Faces',objectHandle.GEOMETRY.faces,...
                    'FaceColor',SIM.OBJECTS(entity).colour,...
                    'EdgeColor',DATA.figureProperties.edgeColor,...
                    'EdgeAlpha',DATA.figureProperties.edgeAlpha,...
                    'FaceLighting',DATA.figureProperties.faceLighting,...
                    'FaceAlpha',DATA.figureProperties.faceAlpha,...
                    'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(entity) = plot3(ax,globalStates(1,frame),globalStates(2,frame),globalStates(3,frame));
                % ASSIGN MARKER PROPERTIES
                markerHandle(entity).Marker = SIM.OBJECTS(entity).symbol;
                markerHandle(entity).MarkerSize = DATA.figureProperties.markerSize;
                markerHandle(entity).MarkerFaceColor = SIM.OBJECTS(entity).colour;
                markerHandle(entity).MarkerEdgeColor = DATA.figureProperties.markerEdgeColor;
                markerHandle(entity).Color = SIM.OBJECTS(entity).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(entity).LineStyle = DATA.figureProperties.lineStyle;
            traceHandle(entity).LineWidth = DATA.figureProperties.lineWidth;
            traceHandle(entity).Color = SIM.OBJECTS(entity).colour;
        end
        % ANNOTATION WITH THE CURRENT TIME
        clockHandle = annotation('textbox',[0.025 0.025 0.15 0.06],'String','Time:','FitBoxToText','off');
        set(clockHandle,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontWeight',DATA.figureProperties.fontWeight);
        % Title
        title(ax,...
            titleString,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.titleFontSize,...
            'FontSmoothing','on');
        % X-Label
        xlabel(ax,'x (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Y-Label
        ylabel(ax,'y (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Z-Label
        zlabel(ax,'z (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Axes
        set(ax,...
            'TickLabelInterpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on',...
            'Color',DATA.figureProperties.axesColor,...
            'GridLineStyle','--',...
            'GridAlpha',0.25,...
            'GridColor','k');

        axisLimit = (DATA.figureProperties.objectMaximalRadii + DATA.figureProperties.maxAbsPosition);
        xlim(ax,[-axisLimit,axisLimit]);
        ylim(ax,[-axisLimit,axisLimit]);
        zlim(ax,[-axisLimit,axisLimit]);
        axis manual;
        % Legend
        legend(markerHandle,DATA.figureProperties.legendEntries,...
            'Location','northeastoutside',...
            'interpreter',DATA.figureProperties.interpreter);
        view([viewPoint inclinationAngle]);     % Set initial view angle
        grid on; box on; hold off; grid minor;
    else
        % CONTINUED FRAMES
        % FOR EACH PLOT ENTITY
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT HANDLE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            globalStates = eval(variableLabel);                                 % Evalutate object states
            % CHECK IF UPDATE IS NECESSARY
            if any(isnan(globalStates(:,frame)))                                % Object is static, freeze its position
                continue
            end
            
            % HANDLE DIFFERENT REPRESENTATIONS
            if numel(objectHandle.GEOMETRY.vertices) > 0
                % GET THE GLOBAL POSE AT STEP
                [R_frame] = OMAS_geometry.quaternionToRotationMatrix(globalStates(7:end,frame));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',objectHandle.GEOMETRY.vertices*R_frame + globalStates(1:3,frame)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',globalStates(1,frame));
                set(markerHandle(entity),'YData',globalStates(2,frame));
                set(markerHandle(entity),'ZData',globalStates(3,frame));
            end
            
            % UPDATE TRACE DATA
            if frame <= DATA.figureProperties.tailLength/SIM.TIME.dt
                % TRANSITIONING PERIOD                    
                tailTrace = globalStates(1:3,1:frame);                     % All points upto the step position
            elseif (frame - DATA.figureProperties.tailLength/SIM.TIME.dt) > 0
                % IF THE TRAIL IS NOW A SUBSET
                tailTrace = globalStates(1:3,(frame-(DATA.figureProperties.tailLength/SIM.TIME.dt)):frame);        % All points upto the step position
            else
                tailTrace = globalStates(1:3,frame);        
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', tailTrace(1,:));
            set(traceHandle(entity),'YData', tailTrace(2,:));
            set(traceHandle(entity),'ZData', tailTrace(3,:));
        end
    end
    % UPDATE TIMESTAMP ANNOTATION
    set(clockHandle,'String',sprintf('Time: %ss',num2str(frameTimeVector(frame))));
    % FORCE IMAGE WRITE
    drawnow();
    
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    im = frame2im(F(frame));
    [imind,cm] = rgb2ind(im,256);
    
    % APPEND THE FRAMES TO GIF
    if frame == 1
        imwrite(imind,cm,filePath,'gif','Loopcount',inf,'DelayTime',(1/DATA.figureProperties.fps));
    else
        imwrite(imind,cm,filePath,'gif','WriteMode','append','DelayTime',(1/DATA.figureProperties.fps));
    end
%     viewPoint = viewPoint + viewRate;
end
% FIGURE COMPLETE
currentFigure = currentFigure + 1;
% CLEAN UP
clearvars -except currentFigure figureHandle
end
% GET THE 3D TRAJECTORY TRAIS AS A VIDEO [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricAvi(SIM,objectIndex,DATA,currentFigure)
% This function generates a new 3D trajectory figure be exporting the
% trajectory data to the output directory before generating the figure.

% CONFIGURE THE PLOT ATTRIBUTES
fileName = strcat(SIM.outputPath,'isometricFigure');
figureHandle = figure('Name','OpenMAS isometric timelapse (AVI)');
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));

% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
if DATA.figureProperties.publish
    set(figureHandle,'MenuBar','none');
    set(figureHandle,'ToolBar','none');
    % set(figureHandle,'Visible','off');
end
% CREATE THE AXES IN THE FIGURE
ax = axes(figureHandle);
% ASSUME THE objectIndex is in the same order as the SIM.OBJECTS

% GET THE SIZE OF THE PLOTTED SET
% The data has already been pre-formatted so that the number of states
% aligns with the number of expected frames in the annimation.
% IMPORT THE OBJECT DATA FROM TEMP FILE
load([SIM.outputPath,SIM.systemFile]);
framesToPlot = size(eval(sprintf('objectID%d',SIM.OBJECTS(1).objectID)),2);   % The variable name in the workspace

% PREPARE THE 'AVI' GENERATOR
vidFile = VideoWriter(fileName);
vidFile.FrameRate = DATA.figureProperties.fps;                           % Set the videos fps to match the sample rate
vidFile.Quality = 50;
open(vidFile);

% ax.NextPlot = 'replaceChildren';
inclinationAngle = 62;
viewPoint = -10; 
% viewRate = 0.05;
for frame = 1:framesToPlot
    % Must be done on a step by step base to get the annimations correct.
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % UPDATE OBJECT HANDLES
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT FOR REFERENCE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            objectStates = eval(variableLabel);                              % Evalutate object states        
            % TRACE INITIAL POSITIONS
            traceHandle(entity) = plot3(ax,objectStates(1,frame),objectStates(2,frame),objectStates(3,frame));
            % THE 
            if numel(objectHandle.GEOMETRY.vertices) > 0
               % GET THE GLOBAL POSE AT STEP
               [R_frame] = OMAS_geometry.quaternionToRotationMatrix(objectStates(7:end,frame));
               % BUILD MARKER
               markerHandle(entity) = patch(ax,...
                   'Vertices',objectHandle.GEOMETRY.vertices*R_frame + objectStates(1:3,frame)',...
                   'Faces',objectHandle.GEOMETRY.faces,...
                   'FaceColor',SIM.OBJECTS(entity).colour,...
                   'EdgeColor',DATA.figureProperties.edgeColor,...
                   'EdgeAlpha',DATA.figureProperties.edgeAlpha,...  
                   'FaceLighting',DATA.figureProperties.faceLighting,...
                   'FaceAlpha',DATA.figureProperties.faceAlpha,...
                   'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(entity) = plot3(ax,objectStates(1,frame),objectStates(2,frame),objectStates(3,frame));
                % ASSIGN MARKER PROPERTIES
                markerHandle(entity).Marker = SIM.OBJECTS(entity).symbol;
                markerHandle(entity).MarkerSize = DATA.figureProperties.markerSize;
                markerHandle(entity).MarkerFaceColor = SIM.OBJECTS(entity).colour;
                markerHandle(entity).MarkerEdgeColor = DATA.figureProperties.markerEdgeColor;
                markerHandle(entity).Color = SIM.OBJECTS(entity).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(entity).LineStyle = DATA.figureProperties.lineStyle;
            traceHandle(entity).LineWidth = DATA.figureProperties.lineWidth;
            traceHandle(entity).Color = SIM.OBJECTS(entity).colour;
        end
        % ANNOTATION WITH THE CURRENT TIME
        clockHandle = annotation('textbox',[0.025 0.025 0.15 0.06],'String','Time:','FitBoxToText','off');
        set(clockHandle,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontWeight',DATA.figureProperties.fontWeight);
        % Title
        title(ax,...
            titleString,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.titleFontSize,...
            'FontSmoothing','on');
        % X-Label
        xlabel(ax,'x (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Y-Label
        ylabel(ax,'y (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Z-Label
        zlabel(ax,'z (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Axes
        set(ax,...
            'TickLabelInterpreter',DATA.figureProperties.interpreter,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on',...
            'Color',DATA.figureProperties.axesColor,...
            'GridLineStyle','--',...
            'GridAlpha',0.25,...
            'GridColor','k');
%         xlim(ax,[DATA.figureProperties.axisMinimums(1),DATA.figureProperties.axisMaximums(1)]);
%         ylim(ax,[DATA.figureProperties.axisMinimums(2),DATA.figureProperties.axisMaximums(2)]);
%         zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);
        axisLimit = (DATA.figureProperties.objectMaximalRadii + DATA.figureProperties.maxAbsPosition);
        xlim(ax,[-axisLimit,axisLimit]);
        ylim(ax,[-axisLimit,axisLimit]);
        zlim(ax,[-axisLimit,axisLimit]);
        axis manual;
        % Legend
        legend(markerHandle,...
            DATA.figureProperties.legendEntries,...
            'Location','northeastoutside',...
            'Interpreter',DATA.figureProperties.interpreter);
%         set(ax,'OuterPosition', [.1, .2, .9, .6]); % [xLeft, yBottom, width, height]
        view([viewPoint inclinationAngle]);     % Set initial view angle
        grid on; grid minor;
        box on; hold off; 
    else
        % CONTINUED FRAMES
        % FOR EACH PLOT ENTITY
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT FOR REFERENCE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            objectStates = eval(variableLabel);                                 % Evalutate object states
            % CHECK IF UPDATE IS NECESSARY
            if any(isnan(objectStates(:,frame)))                                % Object is static, freeze its position
                continue
            end
            % HANDLE DIFFERENT REPRESENTATION
            if numel(objectHandle.GEOMETRY.vertices) > 0
                % GET THE GLOBAL POSE AT STEP
                [R_frame] = OMAS_geometry.quaternionToRotationMatrix(objectStates(7:end,frame));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',objectHandle.GEOMETRY.vertices*R_frame + objectStates(1:3,frame)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',objectStates(1,frame));
                set(markerHandle(entity),'YData',objectStates(2,frame));
                set(markerHandle(entity),'ZData',objectStates(3,frame));
            end
            % UPDATE TRACE DATA
            if frame <= DATA.figureProperties.tailLength/SIM.TIME.dt
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(1:3,1:frame);                     % All points upto the step position
            elseif (frame - DATA.figureProperties.tailLength/SIM.TIME.dt) > 0
                % IF THE TRAIL IS NOW A SUBSET
                tailTrace = objectStates(1:3,(frame-(DATA.figureProperties.tailLength/SIM.TIME.dt)):frame);        % All points upto the step position
            else
                tailTrace = objectStates(1:3,frame);        
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', tailTrace(1,:));
            set(traceHandle(entity),'YData', tailTrace(2,:));
            set(traceHandle(entity),'ZData', tailTrace(3,:));
        end
    end
    % UPDATE TIMESTAMP ANNOTATION
    set(clockHandle,'String',sprintf('Time: %ss',num2str(frameTimeVector(frame))));
    % FORCE IMAGE WRITE
    drawnow;
    % COLLECT FRAMES    
    writeVideo(vidFile,getframe(figureHandle));       % Write frame to video
    % INCREMENT
%     viewPoint = viewPoint + viewRate;
end
close(vidFile);
currentFigure = currentFigure + 1;
% CLEAN UP
clearvars -except currentFigure figureHandle
end

%% /////////////////////////// TIMING FIGURES /////////////////////////////
% COMPUTATION TIME FIGURES
function [currentFigure,figureHandle] = get_computationTimes(SIM,objectIndex,DATA,currentFigure)
% This function gets the computation time figures for all agents using the
% objectIndex.DATA fields.
% INPUTS:
% SIM     - Local copy of the META structure
% OBJECTS - The object class objects.
% currentFigure - The current figure number
% OUTPUTS:
% currentFigure - Updated figure number
% figureHandle - Handle to the created figure

figureHandle = [];

% SANITY CHECK 1
if DATA.totalAgents == 0
    warning('No agent data available.'); 
    return 
else
    % LOGICALLY SELECT THE AGENTS FROM THE OBJECT INDEX 
    agentObjectIndex = objectIndex([SIM.OBJECTS.type] == OMAS_objectType.agent); % The agents themselves  
end

iter = 0;
for ID1 = 1:SIM.totalAgents
    % TRY TO GET THE REQUIRED PROPERTIES
    try
       agentObjectIndex{ID1}.DATA.dt;
    catch
       iter = iter + 1; 
    end
end

% SANITY CHECK 2
if iter == SIM.totalAgents
    warning('No agents with defined algorithm timeseries in "DATA.dt".');
    return 
end

% CONFIGURE THE PLOT ATTRIBUTES
figurePath = strcat(SIM.outputPath,'computationTimes');
figureHandle = figure('Name','OpenMAS computation timeseries');
ax = axes(figureHandle);
set(figureHandle,'Position',DATA.figureProperties.windowSettings);         % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
hold on;
for ind = 1:SIM.totalAgents
    % GET OBJECT DATA
    entity = objectIndex{ind};
    % GET THE EQUIVALENT SIM OBJECT
    SIMobject = SIM.OBJECTS([SIM.OBJECTS.objectID] == entity.objectID);
    % GENERATE FIGURES
    plot(ax,DATA.timeVector(1:SIM.TIME.endStep),entity.DATA.dt*1E3,...
         'Color',SIMobject.colour,...
         'LineWidth',DATA.figureProperties.lineWidth,...
         'DisplayName',sprintf('[ID-%s] %s',num2str(SIMobject.objectID),SIMobject.name));
end

% Title
title(ax,...
    sprintf('Agent Computation times over a period of %ss',num2str(SIM.TIME.endTime)),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.titleFontSize,...
    'FontSmoothing','on');
% X-Label
xlabel(ax,'t (s)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Y-Label
ylabel(ax,'Computation Time (ms)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on');
% Legend
legend(ax,...
    'Location','northeastoutside',...
    'Intepreter',DATA.figureProperties.interpreter);
% Axes
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor,...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
grid on; box on; grid minor;
% Show the timestep difference in the figure
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):SIM.TIME.dt:ax.XAxis.Limits(2);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);      
% PUBLISH TO PDF
if DATA.figureProperties.publish
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0');
end

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;   
end

% NOTES:
% WAIT FOR INPUT BEFORE CLOSING ALL
%         if ~SIM.monteCarloMode
%             input(sprintf('\n[%s]\tPress enter to clear tray and exit.\n',SIM.phase));
%         end