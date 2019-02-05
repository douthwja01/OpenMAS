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