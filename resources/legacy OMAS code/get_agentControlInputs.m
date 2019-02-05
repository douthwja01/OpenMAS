% AGENT CONTROL INPUT TRAJECTORIES
function [currentFigure,figureHandle] = get_agentControlInputs(SIM,objectIndex,DATA,currentFigure)
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
agentObjectIndex = objectIndex([SIM.OBJECTS.type] == OMAS_objectType.agent); % The agents themselves

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