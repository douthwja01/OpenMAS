% GET THE GLOBAL TRAJECTORY DATA FOR AN INDIVIDUAL OBJECT
function [currentFigure,figureHandle] = GetFigure_ObjectTrajectory(SIM,DATA,currentFigure,objectNum)

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