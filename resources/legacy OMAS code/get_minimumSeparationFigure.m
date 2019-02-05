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
figurePath = strcat(SIM.outputPath,'minimumSeparations');    

% FIGURE META PROPERTIES
figureHandle = figure('Name','OpenMAS point of closest approach'); 
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.backGroundColor);           % Background colour 

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
maxDisplayObjects = 10;

% MAX-MATRIX > MAXIMUM SEPERATIONS [IDA by IDB]
maxABMatrix = max(separationTimeSeries,[],3);
maxABMatrix(isinf(maxABMatrix)) = NaN;

[h_maxVals,~] = max(maxABMatrix,[],2);       % The minimum seperations, horezontal indicies
maxProximity = max(h_maxVals);               % Maximal seperation

% MIN-MATRIX > MINIMUM SEPERATIONS [IDA by IDB]
minABAxes = collidableIDs;
minABMatrix = min(separationTimeSeries,[],3);                          % Minimum seperations over the timeseries

% minABMatrix is a matrix of closest interactions between all objects A
% vs B.
[h_val,h_ind] = min(minABMatrix,[],2); % The minimum seperations, horezontal indicies

% DETERMINE UNIQUE INTERACTIONS
[ord,ord_ind] = sort(h_val);           % Sort minimums by magnitude
[~,unique_ind,~] = unique(ord,'rows');

uniqueIndOrder = ord_ind(unique_ind);  % Unique, ordered proximities

% THE ID's OF THE CLOSEST (UNIQUE) OBJECTS
closestIDpairs = [minABAxes',h_ind];
orderedClosestPairs = closestIDpairs(uniqueIndOrder,:);
% GET THE LIMITED NUMBER OF ID PAIRS TO PLOT
if size(orderedClosestPairs,1) > maxDisplayObjects
    plotableInteractions = orderedClosestPairs(1:maxDisplayObjects,:);
else
    plotableInteractions = orderedClosestPairs;
end

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
legendEntries{legendCounter} = 'Collision Boundary';

% ADD FINAL PLOT ATTRIBUTES
titleString = sprintf('Separation of the closest %s Objects',num2str(maxDisplayObjects));
title(ax,titleString,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.titleFontSize,'fontSmoothing','On');
xlabel(ax,'t(s)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'fontSmoothing','On');
ylabel(ax,'Separation(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'fontSmoothing','On');
set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'fontSmoothing','On');
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
ylim(ax,[0 (0.5*maxProximity)]);
% legHandle = legend(ax,legendEntries,'Location','northeastoutside');
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