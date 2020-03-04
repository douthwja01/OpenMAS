% AGENT CONTROL INPUT TRAJECTORIES
function [currentFigure,figureHandle] = GetFigure_inputs(SIM,objectIndex,DATA,currentFigure)
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
    warning('No agents have input timeseries stored in (agent).DATA.inputs. Use the "DATA.inputNames" and "DATA.inputs" fields if desired.');
    return 
end

% Declare title string for figure  
figurePath = strcat(SIM.outputPath,'inputs');

% FIGURE META PROPERTIES
figureHandle = figure('Name','OpenMAS control inputs');
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.85]);
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
    inputTrajectories = evalAgent.DATA.inputs;
    inputCount = size(inputTrajectories,1);
    
    % BUILD THE SUB-PLOT
    plotCellB = double(ID1*plotCellWidth);                                         % The end of the plot
    plotLocation = subplot(double(DATA.totalAgents),plotCellWidth,[plotCellA plotCellB]);
    set(plotLocation,...
        'Color',DATA.figureProperties.axesColor,...
        'GridLineStyle','--',...
        'GridAlpha',0.25,...
        'GridColor','k');
    hold on;
    
    % LEGEND LABELS
    legendEntries = cell(inputCount,1);   
    numericLabelling = 1;
    if isfield(evalAgent.DATA,'inputNames') 
        if length(evalAgent.DATA.inputNames) ~= inputCount
            warning('Incorrect number of input labels for agent %s, reverting to numeric labelling',simName);
            numericLabelling = 1;
        else
            % FETCH THE INPUT LABEL
            for entry = 1:size(inputTrajectories,1)
                legendEntries{entry} = evalAgent.DATA.inputNames{entry};
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
        'fontname',DATA.figureProperties.fontName,...
        'fontweight',DATA.figureProperties.fontWeight,...
        'fontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','On');
    % Axes
    set(plotLocation,...
        'TickLabelInterpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
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
        'fontname',DATA.figureProperties.fontName,...
        'Interpreter',DATA.figureProperties.interpreter);
    
    grid on; box on; grid minor;
    % Show the timestep difference in the figure
    %plotLocation.XAxis.MinorTickValues = plotLocation.XAxis.Limits(1):SIM.TIME.dt:plotLocation.XAxis.Limits(2);
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if plotCellA == 1
        % Title
        title(plotLocation,sprintf('Agent control trajectories over a period of %ss',num2str(SIM.TIME.endTime)),...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontsize',DATA.figureProperties.titleFontSize);               % Append title to first subplot
    end
    % Prevent overlap of x-label
    if ID1 == DATA.totalAgents
        % X-Label
        xlabel(plotLocation,'t (s)',...
            'Interpreter', DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','On');                                         
    else
        set(plotLocation,'XTickLabel',[]); 
    end
    % Define the x-limits
    xlim([SIM.TIME.startTime,SIM.TIME.endTime]);                           % Ensures plot alignment/sizing
    plotCellA = plotCellA + plotCellWidth; % Move to next subplot location  
end
hold off;

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