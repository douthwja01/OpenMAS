%% GET THE AGENT SEPERATION SEPERATION TIMESERIES FIGURE
function [currentFigure,figureSet] = GetFigure_avoidanceRates(SIM,objectIndex,DATA,currentFigure)
% Draws the separations between the 'subject object', its associated waypoints 
% and obstacles.

figureSet = [];
if SIM.totalObjects < 2 || SIM.totalAgents == 0
    warning('There must be at least two agents that are collidable to plot avoidance data.\n');
    return
end

% Assemble the seperation figures
agentIDs = [SIM.OBJECTS([SIM.OBJECTS.type] == OMAS_objectType.agent).objectID];
n = 0;
for i = 1:SIM.totalAgents
    agentObject = objectIndex{SIM.globalIDvector == agentIDs(i)};           % Get the agent in question
    
    if isempty(agentObject.DATA) % Check the agent has input trajectory data (may not be recorded)
        warning('Agent %s has no input trajectory data stored in (agent)."DATA.input".',agentObject.name);
        continue;
    else
        n = n + 1;
    end 
    
    % Generate the figure
    [currentFigure,figureSet(n)] = BuildAvoidanceFigure_local(...
        SIM,...
        DATA,currentFigure,...
        agentObject);
end

% If there are no figures
if n < 1
   return 
end

% ASSEMBLE TABBED FIGURE
windowHandle = GetTabbedFigure(figureSet,'OpenMAS avoidance overview');
set(windowHandle,'Position', DATA.figureProperties.windowSettings);         % Maximise the figure in the tab
savefig(windowHandle,[SIM.outputPath,'avoidance-overview']);                % Save the output figure
end

%% PLOT FROM LOCAL STATES
function [currentFigure,figureHandle] = BuildAvoidanceFigure_local(SIM,DATA,currentFigure,subjectOBJ)

% Get the subject data
subjectMETA   = SIM.OBJECTS(SIM.globalIDvector == subjectOBJ.objectID);
subjectStates = OMAS_getTrajectoryData_mex(...
    DATA.globalTrajectories,...
    SIM.globalIDvector,...
    subjectOBJ.objectID,inf);

% Second object references
collidableMETA  = SIM.OBJECTS([SIM.OBJECTS.hitBox] ~= OMAS_hitBoxType.none);         % Only collidable objects
collidableMETA  = collidableMETA([collidableMETA.type] ~= OMAS_objectType.waypoint); % That are not waypoints
collidableMETA  = collidableMETA([collidableMETA.objectID] ~= subjectMETA.objectID); % Not the same object 

% Figure properties
figureHandle = figure('Name',sprintf('Avoidance trajectories of %s',subjectMETA.name));
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.82]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
plotCellA = 1; plotCellWidth = 4;                                          % The width of each figure, the start of the plot

maxSeparation = 0;  % For x-axis scaling
for n = 1:length(collidableMETA)
    % Get the obstacle data
    obstacleMETA   = collidableMETA(n);
    obstacleStates = OMAS_getTrajectoryData_mex(...
        DATA.globalTrajectories,...
        SIM.globalIDvector,...
        obstacleMETA.objectID,...
        inf);
    % Separation constraint
    separationBoundary = obstacleMETA.radius + subjectMETA.radius;
    
    % Calculate the separation timeseries
    timeseries_separation = NaN(1,SIM.TIME.endStep);
    timeseries_speed      = NaN(1,SIM.TIME.endStep);
    timeseries_heading    = NaN(3,SIM.TIME.endStep);
    eta0 = OMAS_geometry.quaternionToEulers(subjectStates(7:10,1));        % Initial heading
    for step = 1:SIM.TIME.endStep
        timeseries_separation(1,step)   = norm(obstacleStates(1:3,step) - subjectStates(1:3,step));
        timeseries_speed(1,step)        = norm(subjectStates(4:6,step));
        timeseries_heading(:,step)      = OMAS_geometry.quaternionToEulers(subjectStates(7:10,step)) - eta0;
    end
    timeseries_heading = timeseries_heading*180/pi;
    
    % Calculate the maximum separation (for x-axis scaling)
    d_max = max(timeseries_separation,[],2);                               
    if d_max > maxSeparation
        maxSeparation = d_max;
    end
    % Arrange for plotting
    timeseries_avoidance = vertcat(...
        timeseries_speed,...
        timeseries_heading(2,:),...
        timeseries_heading(3,:));                                          % Can neglect roll (irrelevant for avoidance) 
    
    %% Begin building the subplot
    plotCellB = double(n*plotCellWidth);                                    % The end of the plot
    plotLocation = subplot(length(collidableMETA),plotCellWidth,[plotCellA plotCellB]);
    set(plotLocation,...
        'GridLineStyle','--',...
        'GridAlpha',0.25,...
        'GridColor','k',...
        'box','on',...
        'XGrid','on',...
        'YGrid','on');
    set(plotLocation,'XDir','reverse');
    
    set(plotLocation,...
    'TickLabelInterpreter',DATA.figureProperties.interpreter,...
    'FontName',DATA.figureProperties.fontName,...
    'FontWeight',DATA.figureProperties.fontWeight,...
    'FontSize',DATA.figureProperties.axisFontSize,...
    'FontSmoothing','on',...
    'Color',DATA.figureProperties.axesColor);

    hold on;
    
    % PLOT THE SEPERATION DATA ON CURRENT FIGURE
    plot(plotLocation,...
         timeseries_separation,...
         timeseries_avoidance,...
        'LineStyle','-',...
        'LineWidth',DATA.figureProperties.lineWidth); 
    set(plotLocation,'Xlim',[0 10]);
    set(plotLocation,'Xdir','reverse'); 
%    set(plotLocation,'XTick',...
%        plotLocation.XAxis.Limits(1):1:plotLocation.XAxis.Limits(2));
    
    % Y-Label
    ylabel(plotLocation,...
        sprintf('w.r.t [ID-%s]',num2str(obstacleMETA.objectID)),...
        'Interpreter',DATA.figureProperties.interpreter,...
        'fontname',DATA.figureProperties.fontName,...
        'fontweight',DATA.figureProperties.fontWeight,...
        'fontSize',DATA.figureProperties.axisFontSize,...
        'FontSmoothing','On');
    
    % Add separation constraint line
    l = xline(plotLocation,separationBoundary);
    set(l,...
        'Color','b',...
        'LineWidth',DATA.figureProperties.lineWidth,...
        'LineStyle','-.',...
        'Label',"Collision Boundary");

    % Legend
    legend(plotLocation,...
        {'$v (m/s)$';'$\theta (deg)$';'$\psi (deg)$'},...
        'Location','eastoutside',...
        'fontname',DATA.figureProperties.fontName,...
        'FontSize',DATA.figureProperties.axisFontSize,...
        'Interpreter',DATA.figureProperties.interpreter);
    
    % Titles
    if n == 1                       % Add title on first subplot
        titleStr = sprintf("Agent %s's avoidance trajectories over the %ss period.",...
            subjectMETA.name,num2str(SIM.TIME.endTime));
        title(plotLocation,titleStr,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontsize',DATA.figureProperties.titleFontSize);                % Append title to first subplot
    end
    if n == length(collidableMETA)  % Prevent overlap of x-label
        xlabel(plotLocation,'Separation (m)',...                            % X-Label
            'Interpreter', DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'fontweight',DATA.figureProperties.fontWeight,...
            'fontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','On');                                         
    else
        set(plotLocation,'XTickLabel',[]); 
    end  
    plotCellA = plotCellA + plotCellWidth; % Move to next subplot location  
end

% General axes settings
set(findall(figureHandle,'-property','XLim'),'XLim',[0 maxSeparation]);
% set(findall(figureHandle,'-property','Units'),'Units','Inches');

figurePath = strcat(...
    SIM.outputPath,...
    sprintf('avoidance-parameters-id-%d-%s',subjectMETA.objectID,subjectMETA.name)); 
% SAVE THE OUTPUT FIGURE
savefig(figureHandle,figurePath);   

% Publish as .pdf if requested
if DATA.figureProperties.publish
	GetFigurePDF(figureHandle,figurePath);   
end

currentFigure = currentFigure + 1;
end