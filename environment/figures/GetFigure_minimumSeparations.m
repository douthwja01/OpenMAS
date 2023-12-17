% GET THE CLOSEST PROXIMITY FIGURES
function [currentFigure,figureHandle] = GetFigure_minimumSeparations(SIM,DATA,currentFigure)
% This function simply generates plots for a given number of the nearest
% misses / smallest separations.

% Defaults
maximumToDisplay = 10;
figureHandle = []; 

% Get the total number of collidable objects
collidableMETA = SIM.OBJECTS([SIM.OBJECTS.hitBox] ~= OMAS_hitBoxType.none);         % Must be collidable
collidableMETA = collidableMETA([collidableMETA.type] ~= OMAS_objectType.waypoint); % Must not be way-points

% Input sanity check #1
if numel(collidableMETA) < 2
    warning('There must be at least two objects that are both collidable and not way-points.');
    return
end

%% We need to cross-reference each object
ind = 1; conflictStruct = struct();
for i = 1:numel(collidableMETA)
    % Retrieve the states of object A
%     states_A = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,collidableMETA(i).objectID,inf);
    states_A = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,collidableMETA(i).objectID,inf);
        
    states_A = states_A(:,1:SIM.TIME.endStep);
    for j = (i+1):numel(collidableMETA)       
        % Get the state trace for object B
%         states_B = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,collidableMETA(j).objectID,inf);
        states_B = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,collidableMETA(j).objectID,inf);
                
        states_B = states_B(:,1:SIM.TIME.endStep);
        % Calculate the separation across the time period
        centroidSeparations = states_B(1:3,:) - states_A(1:3,:);
        centroidSeparations = sqrt(sum(centroidSeparations.^2,1));
        % Retain the structure as a vector
        conflictStruct(ind).objectID_A = collidableMETA(i).objectID;
        conflictStruct(ind).objectID_B = collidableMETA(j).objectID;
        conflictStruct(ind).mag = centroidSeparations;
        conflictStruct(ind).min = min(centroidSeparations);
        conflictStruct(ind).max = max(centroidSeparations);
        % Calculate the critical separation
        conflictStruct(ind).constraint = collidableMETA(i).radius + collidableMETA(j).radius;
        % Increment the reference vector
        ind = ind + 1;
    end
end
% Reorder on minimal proximity
[~,ind] = sort([conflictStruct.min]);
if numel(ind) > maximumToDisplay
    ind = ind(1:maximumToDisplay);
end
% Redefine the conflict array
conflictStruct = conflictStruct(ind);

%% Define the figure
figurePath = strcat(SIM.outputPath,'minimum-separations');   
figureHandle = figure('Name','OpenMAS point of closest approach'); 
setappdata(figureHandle,'SubplotDefaultAxesLocation', [0.1, 0.1, 0.80, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
ax = axes(figureHandle);                                                   % Begin generating the figure                                   
grid on; box on; hold on; grid minor;

% Generate the line references
legendEntries = cell(numel(conflictStruct),1);
for i = 1:numel(conflictStruct)
    % The associate META objects
    META_A = collidableMETA([collidableMETA.objectID] == conflictStruct(i).objectID_A);  
    META_B = collidableMETA([collidableMETA.objectID] == conflictStruct(i).objectID_B);
    % Create the legend entries
    legendEntries{i} = sprintf('%s[ID-%d] w.r.t. %s[ID-%d]',META_A.name,META_A.objectID,META_B.name,META_B.objectID);
    % Generate the plot
    l = plot(ax,DATA.timeVector(:,1:SIM.TIME.endStep),conflictStruct(i).mag(:,1:SIM.TIME.endStep));
    set(l,'Color',META_A.colour);
    set(l,'LineStyle','-');
    set(l,'LineWidth',DATA.figureProperties.lineWidth);
end

% [TO-DO]
% - Handle legend entries for multiple collision constraints

try
    % Generate the reference line for the last conflict
    c = refline(ax,0,conflictStruct(i).constraint);                            % Adds a reference line with slope m and intercept b to the current axes.
    set(c,'Color','k');
    set(c,'LineStyle','--');
    set(c,'LineWidth',DATA.figureProperties.lineWidth/2);
    % Add legend entries
    legendEntries{i+1} = 'Collision Boundary';
catch ex
    warning("Failed to plot collision condition ref-line reason:\n %s",ex.message);
end

%% Add additional figure refinements 
% Title
title(ax,sprintf('Global minimum separations over a period of %ss',num2str(SIM.TIME.endTime)),...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'fontsize',DATA.figureProperties.titleFontSize,...
    'fontweight',DATA.figureProperties.fontWeight,...
    'fontsmoothing','On');
% X-axes
xlabel(ax,'t (s)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSmoothing','On');
xlim(ax,[SIM.TIME.startTime SIM.TIME.endTime]);
% Y-axes
ylabel(ax,'Separation (m)',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'Fontweight',DATA.figureProperties.fontWeight,...
    'FontSmoothing','On');
ylim(ax,[0 1.1*max([conflictStruct.max])]);                                % Use the maximum recorded separation 
% Add legend entries
legend(ax,legendEntries,...
    'Location','northEastOutside',...
    'Interpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName);
% Axes properties
set(ax,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'fontname',DATA.figureProperties.fontName,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontWeight',DATA.figureProperties.fontWeight,...
    'Color',DATA.figureProperties.axesColor,...
    'FontSmoothing','On',...
    'GridLineStyle','--',...
    'GridAlpha',0.25,...
    'GridColor','k');
% Show the timestep difference in the figure
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):SIM.TIME.dt:ax.XAxis.Limits(2);

% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);   

% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end

% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;    
end
