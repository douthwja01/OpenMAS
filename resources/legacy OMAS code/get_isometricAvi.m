% GET THE 3D TRAJECTORY TRAIS AS A VIDEO
function [currentFigure,figureHandle] = get_isometricAvi(SIM,DATA,currentFigure)
% This function generates a video of the complete simulation trajectories
% over the complete series.

% OUTPUT FILE
fileName = strcat(SIM.outputPath,'isometricVideo','.avi');

% DATA CONTAINERS
xData = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
yData = xData; zData = xData;
xMarker = NaN(1,SIM.totalObjects,SIM.TIME.numSteps); 
yMarker = xMarker; zMarker = xMarker;
% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.endStep
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active
        % DESIRED TRACE LENGTH
        tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length
        if step < activeSteps 
            % IF THE AGENT IS ACTIVE
            markerPosition = objectStates(1:3,step);                       % Step position
            if step > tailStepLength
                % ONLY THE STATES FOR THE TRACE DURATION
                traceEnd = step - (tailStepLength-1);
                tailTrace = objectStates(1:3,traceEnd:step);               % States between the step and (step - tailStepLength)
            else
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(1:3,1:step);                      % All points upto the step position
            end
        else
            % IF THE AGENT IS INACTIVE
            markerPosition = objectStates(1:3,end);                        % Step position
            tailTrace = objectStates(1:3,((activeSteps + 1) - tailStepLength):end);
        end
        
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        xMarker(1,indexValue,step) = markerPosition(1);
        yMarker(1,indexValue,step) = markerPosition(2);
        zMarker(1,indexValue,step) = markerPosition(3); % Coordinate by object by step
        % BUILD MATICES OF TAIL COORDINATES
        traceData = tailTrace';
        xData(1:size(traceData,1),indexValue,step) = traceData(:,1); % Full X coordinate set
        yData(1:size(traceData,1),indexValue,step) = traceData(:,2); % Full X coordinate set        
        zData(1:size(traceData,1),indexValue,step) = traceData(:,3); % Full X coordinate set
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Video Timelapse');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', [DATA.figureProperties.alignment 50 DATA.figureProperties.size (DATA.figureProperties.size-130)]); % [x y width height]
set(figureHandle,'MenuBar','none');
set(figureHandle,'ToolBar','none');
set(figureHandle,'Visible','off');

% CONFIGURE VIDEO SETTINGS
fps = 50; % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.simTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end
% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(SIM.TIME.endStep/numFrames);                         % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end
vidFile = VideoWriter(fileName);
vidFile.FrameRate = fps;                                                   % Set the videos fps to match the sample rate
open(vidFile);
% FRAME SETTINGS
frameSet = moviein(numFrames);                                             % Pre-allocate video object for frames
F(numFrames) = struct('cdata',[],'colormap',[]);
% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(DATA.timeVector(end)));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1; viewPoint = -45;
for step = 1:SIM.TIME.endStep
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % INITIAL TRACE PLOT
        traceHandle = plot3(ax,xData(:,:,step),yData(:,:,step),zData(:,:,step));
        for ind = 1:SIM.totalObjects
            % INITIAL MARKER PLOT
            markerHandle(ind) = plot3(ax,xMarker(1,ind,step),yMarker(1,ind,step),zMarker(1,ind,step));
            % ASSIGN MARKER PROPERTIES
            markerHandle(ind).Marker = SIM.OBJECTS(ind).symbol;
            markerHandle(ind).MarkerSize = DATA.figureProperties.MarkerSize;%*SIM.OBJECTS(ind).radius;
            markerHandle(ind).MarkerFaceColor = SIM.OBJECTS(ind).colour;
            markerHandle(ind).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
            markerHandle(ind).Color = SIM.OBJECTS(ind).colour;
            % ASSIGN TRACE PROPERTIES
            traceHandle(ind).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(ind).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(ind).Color = SIM.OBJECTS(ind).colour;
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('y_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        zlabel('z_{m}','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
                
        xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        grid on; grid minor; box on; hold off; 

        legend(markerHandle,DATA.figureProperties.legendEntries,'location','northeastoutside');
        view([viewPoint 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            % UPDATE MARKER DATA
            set(markerHandle(entity),'XData',xMarker(:,entity,step));
            set(markerHandle(entity),'YData',yMarker(:,entity,step));
            set(markerHandle(entity),'ZData',zMarker(:,entity,step));
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity,step));
            set(traceHandle(entity),'YData', yData(:,entity,step));
            set(traceHandle(entity),'ZData', zData(:,entity,step));
        end
        view([viewPoint 50]);
        viewPoint = viewPoint + 0.05;
    end
    drawnow();
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    writeVideo(vidFile,F(frame));       % Write frame to video
    frame = frame + 1;                  % Move to next frame
end
close(vidFile);
% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end
