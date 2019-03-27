% COMPUTATION TIME FIGURES
function [currentFigure,figureHandle] = GetFigure_ComputationTimes(SIM,objectIndex,DATA,currentFigure)
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