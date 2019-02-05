% GET THE OBJECT SEPERATION TIMESERIES FIGURE
function [currentFigure,figureHandle] = get_agentSeparationFigure(SIM,DATA,currentFigure)
% Draws the seperations of all agents from each other, neglecting waypoint
% objects.

if DATA.totalAgents < 2
	warning('There must be at least one agent to observe separations.');
    figureHandle = [];
    return
end

% Declare title string for figure  
titlestr = sprintf('Agent separation over a period of %ss',num2str(SIM.TIME.endTime));

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Separation Overview');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 0.83]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.backGroundColor);           % Background colour 
plotCellWidth = 4; plotCellA = 1;                                          % The width of each figure, the start of the plot

% FOR EACH AGENT'S PERSPECTIVE
for ID1 = 1:DATA.totalAgents
    % CHECK IF EVALUATION OBJECT IS A WAYPOINT; OMIT IF NECESSARY
    if SIM.OBJECTS(ID1).type == OMAS_objectType.waypoint
        continue
    end
    
    % EXTRACT TIME-STATE DATA FROM THE TRAJECTORY MATRIX
    [objectID1States] = OMAS_getTrajectoryData(DATA,ID1);              % Trajectory data per objectID
    objectID1 = SIM.OBJECTS(ID1).objectID;
    objectID1Name = SIM.OBJECTS(ID1).name;
    ystring = sprintf('%s [ID:%s] (m)',objectID1Name,num2str(objectID1));
    
    % CALCULATE THE ABSOLUTE SEPERATION FROM ALL OTHER OBJECTS OVER TIME
    legendEntries = cell(1,(DATA.totalAgents-1));
    legendCounter = 1; maxSeriesSeperation = 1;                            % Initial y-axis limit
    plotCellB = ID1*plotCellWidth;                                         % The end of the plot
    plotLocation = subplot(DATA.totalAgents,plotCellWidth,[plotCellA plotCellB]);
    hold on;
    
    % MOVE THROUGH OTHER OBJECTS
    for ID2 = 1:DATA.totalObjects
        % DEFINE CONDITIONS FOR PLOTTING
        plotCondition = SIM.OBJECTS(ID1).objectID ~= SIM.OBJECTS(ID2).objectID && ~strcmpi(SIM.OBJECTS(ID2).type,'WAYPOINT');
        if plotCondition
            objectID2 = SIM.OBJECTS(ID2).objectID;                    % Assemble the object legend information
            objectID2Name = num2str(SIM.OBJECTS(ID2).name);
            % Generate legend entry
            legendEntries(legendCounter) = {sprintf('[ID:%s] %s',num2str(objectID2),objectID2Name)};
            % Container for the timeseries data
            ABSseperationTimeSeries = zeros(1,size(objectID1States,2));
            % Get comparative state data
            [objectID2States] = OMAS_getTrajectoryData(DATA,ID2);
            
            % Reset collision instance to 0
            collisionInstance = 0;
            for simStep = 1:size(ABSseperationTimeSeries,2)
                % CALCULATE THE ABSOLUTE OBSTACLE SEPERATION TIMESERIES
                collisionCondition = SIM.OBJECTS(ID1).radius + SIM.OBJECTS(ID2).radius;                         % Define the collision condition (physical seperation)                      
                objectSeperation = sqrt(sum((objectID1States(1:3,simStep) - objectID2States(1:3,simStep)).^2)); % Define object absolute (positional seperation)  
                objectSeperation = objectSeperation - collisionCondition;                                       % Subtract the physical seperation requirement from the positional seperation
                
                if 0 >= objectSeperation && collisionInstance == 0
                    % OBJECT COLLISION INSTANCE
                    collisionInstance = 1;
                    % ADD COLLISION ANNOTATION
                    annoString = strcat('Collision-',objectID2Name);
                    xCoord = DATA.timeVector(simStep);
                    yCoord = objectSeperation;
                    text(plotLocation,(xCoord),(yCoord),'|','fontsize',20,'Color',SIM.OBJECTS(ID2).colour); % '\uparrow'
                    text(plotLocation,xCoord,(yCoord + 1),annoString);
                end
                ABSseperationTimeSeries(simStep) = objectSeperation;
            end
            
            % Calculate seperation axis limit
            maxSeperation = max(ABSseperationTimeSeries); % Get the maximum seperation value
            if maxSeriesSeperation < maxSeperation 
                maxSeriesSeperation = maxSeperation;
            end
            
            %% PLOT THE SEPERATION DATA ON CURRENT FIGURE
            plot(plotLocation,DATA.timeVector,ABSseperationTimeSeries,...
                        'LineStyle','-',...
                        'LineWidth',DATA.figureProperties.LineWidth,...                
                        'Color',SIM.OBJECTS(ID2).colour);
            ylabel(ystring,'fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
            ylim([(-0.1*maxSeriesSeperation) (1.1*maxSeriesSeperation)]);
            grid on; box on;
            set(gca,'FontSize',DATA.figureProperties.axisFontSize,...
                  'fontWeight',DATA.figureProperties.fontWeight);
            legendCounter = legendCounter + 1;                             % Increment legend entry
        end
    end
    
    % ADD FIGURE TITLE TO FIRST PLOT HEADER
    if ID1 == 1
        title(titlestr,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);                                           % Append title to first subplot
    end
    % PREVENT X LABEL OVERLAP
    if ID1 == DATA.totalAgents
        xlabel(plotLocation,'t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize,'FontSmoothing','on');
%         xlabel('t_{s}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
    else
        set(plotLocation,'XTickLabel',[]); % 'XTick',[]
    end
    % ADD LEGEND
    legend(legendEntries,'Location','eastoutside');
           
    % ADD ANNOTATION THE FOR THE VISUAL HORIZON
%     line([SIM.TIME.startTime SIM.TIME.endTime],...
%          [SIM.OBJECTS(ID1).detectionRange SIM.OBJECTS(ID1).detectionRange],...
%          'LineStyle','--','color','bl');
%     text((SIM.TIME.startTime + SIM.TIME.dt),(SIM.OBJECTS(ID1).detectionRange - 5),...
%          'Visual horizon','color','bl');   % Annotate line as visual horizon

    plotCellA = plotCellA + plotCellWidth; % Move to next subplot location
end
hold off;
% SAVE THE OUTPUT FIGURE
savefig(figureHandle,strcat(SIM.outputPath,'separations.fig'));            % As matlab figure 
saveas(figureHandle,strcat(SIM.outputPath,'separations'),'epsc');          % As postscript
% INCREMENT THE FIGURE INDEX
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;      
end