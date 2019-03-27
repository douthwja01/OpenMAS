% GET THE CLOSEST PROXIMITY FIGURES
function [currentFigure,figureHandle] = GetFigure_MinimumSeparations(SIM,DATA,currentFigure)
% This function generates a figure with the smallest separations observed
% in the simulation and plots their separation. This plot only plots the
% separation between objects that have a form of hitbo (i.e. ~none).

figureHandle = []; % Default figure container

% Input sanity check #1
if ~any([SIM.OBJECTS(:).hitBox] ~= OMAS_hitBoxType.none) 
    warning('There must be at least two collidable objects to plot seperation data.');
    return
end
% Input sanity check #2
if ~any([SIM.OBJECTS(:).type] ~= OMAS_objectType.waypoint) 
    warning('There are no objects that are not waypoints.');
    return
end

% Get the META structures for all hitbox objects
hitBoxOBJECTS = SIM.OBJECTS([SIM.OBJECTS.hitBox] ~= OMAS_hitBoxType.none);   % The objects with valid hitboxes
hitBoxOBJECTS = hitBoxOBJECTS([hitBoxOBJECTS.type] ~= OMAS_objectType.waypoint);

% Input sanity check #3 
if numel(hitBoxOBJECTS) < 2
    warning('Only one collidable object, no separations can be plotted.');
    return
end

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