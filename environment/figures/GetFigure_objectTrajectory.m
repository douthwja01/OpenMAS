% GET THE GLOBAL TRAJECTORY DATA FOR AN INDIVIDUAL OBJECT
function [currentFigure,figureSet] = GetFigure_objectTrajectory(SIM,DATA,currentFigure)

figureSet = [];
for n = 1:DATA.totalObjects
    [currentFigure,figureSet(n)] = BuildObjectTrajectoryFigure(SIM,DATA,currentFigure,n);
end
% ASSEMBLE TABBED FIGURE
windowHandle = GetTabbedFigure(figureSet,'OpenMAS Trajectory Overview');
set(windowHandle,'Position',DATA.figureProperties.windowSettings);          % Maximise the figure in the tab
savefig(windowHandle,[SIM.outputPath,'global-trajectory-overview']);        % Save the output figure
end

% Generate the trajectory figure for the object
function [currentFigure,figureHandle] = BuildObjectTrajectoryFigure(SIM,DATA,currentFigure,objectNum)

% Declare title string for figure  
titlestr = sprintf('Global trajectory data for %s over a period of %ss',SIM.OBJECTS(objectNum).name,num2str(SIM.TIME.endTime));

objectLabel = sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(objectNum).objectID),SIM.OBJECTS(objectNum).name);

% CONFIGURE THE PLOT ATTRIBUTES
figurePath  = strcat(SIM.outputPath,sprintf('global-trajectories-%s',sprintf('id-%s-%s',num2str(SIM.OBJECTS(objectNum).objectID),SIM.OBJECTS(objectNum).name)));
figureHandle = figure('Name',objectLabel);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
set(figureHandle,'Visible','off');
plotCellWidth = 4; plotCellA = 1;                                          % The width of each figure, the start of the plot
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.80]);  

% EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
%[objectStates] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(objectNum).objectID,inf);
[objectStates] = OMAS_getTrajectoryData(DATA.globalTrajectories,SIM.globalIDvector,SIM.OBJECTS(objectNum).objectID,inf);

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
        'fontname',DATA.figureProperties.fontName,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on');
    % Axes
    set(plotLocation,...
        'TickLabelInterpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'Fontweight',DATA.figureProperties.fontWeight,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','on',...
        'Color',DATA.figureProperties.axesColor);
    grid on; box on; grid minor;
    plotLocation.XAxis.MinorTickValues = plotLocation.XAxis.Limits(1):SIM.TIME.dt:plotLocation.XAxis.Limits(2);
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if n == 1
        % Append title to first subplot
        title(plotLocation,titlestr,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'Fontsize',DATA.figureProperties.titleFontSize,...
            'FontSmoothing','on');   
        
    end
    % Prevent overlap of x-label
    if n == stateVectorLength
        % X-axis
        xlabel(plotLocation,'t (s)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
    else
        set(plotLocation,'XTickLabel',[]); 
    end
    % Define the x-limits
    xlim([SIM.TIME.startTime,SIM.TIME.endTime]);                           % Ensures plot alignment/sizing
    % Move to next subplot location 
    plotCellA = plotCellA + plotCellWidth;
end
hold off; 

% SAVE THE OUTPUT FIGURE
set(figureHandle,'Visible','on');                % Make it visable for saving
savefig(figureHandle,figurePath);                % As matlab figure 
% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end
set(figureHandle,'Visible','off');
currentFigure = currentFigure + 1;  % Iterate the plot count
end
