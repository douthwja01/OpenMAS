% GET THE 3D TRAJECTORY TRAILS AS A GIF  [ UPDATED ]
function [currentFigure,figureHandle] = get_isometricGif(SIM,DATA,currentFigure)
% This function generates an animated .gif representing the object
% trajectories over the complete timeseries.

% OUTPUT FILE
filePath = strcat(SIM.outputPath,'isometricFigure.gif');

% DATA CONTAINERS
% We need to build a matrix of states*tailLength*IDs,step
globalTraces = NaN((DATA.figureProperties.tailLength/SIM.TIME.dt),10,DATA.totalObjects,SIM.TIME.numSteps); % Pre-allocate matrix
globalMarker = NaN(10,SIM.totalObjects,SIM.TIME.numSteps); 
% DESIRED TRACE LENGTH
tailStepLength = DATA.figureProperties.tailLength/SIM.TIME.dt;     % Get the tail step length

% BEGIN GENERATING THE MARKER AND TRACE PLOT MATRICES
for step = 1:SIM.TIME.endStep
    for indexValue = 1:DATA.totalObjects
        % GET THE TRAJECTORY DATA FOR THE AGENT
        [objectStates] = OMAS_getTrajectoryData(DATA,indexValue,'valid'); % All valid states for the object      
        activeSteps = size(objectStates,2);                                % The number of steps where object is active

        % IF SIMULATION IS JUST STARTING
        if step < tailStepLength && step < activeSteps
        	% TRANSITIONING PERIOD                    
        	tailTrace = objectStates(:,1:step);                            % All points upto the step position
        % IF THE TRAIL IS NOW A SUBSET
        elseif step <= activeSteps && (step-tailStepLength) > 0
            tailTrace = objectStates(:,(step-tailStepLength):step);        % All points upto the step position
        else
            tailTrace = objectStates(:,step);        
        end
        markerState = objectStates(:,step);                            % Step position
                      
        % COLUMNS OF MATRICES X,Y,Z ARE THE OBJECTS
        globalMarker(:,indexValue,step) = markerState;      
        
        % BUILD MATICES OF TAIL COORDINATES
        % We need to build a matrix of dimensions
        % [tailLength*states*IDs,step], the resulting plot matrices must be
        % of dimensions [state(1)*tailLength,IDs]
        globalTraces(1:size(tailTrace,2),:,indexValue,step) = tailTrace';
    end
end

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS Isometric Timelapse');
% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
% set(figureHandle,'MenuBar','none');
% set(figureHandle,'ToolBar','none');
% set(figureHandle,'Visible','off');


% CONFIGURE VIDEO SETTINGS
fps = 50; % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.endTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end
% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(SIM.TIME.endStep/numFrames);                         % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end

% F(numFrames) = struct('cdata',[],'colormap',[]);
F(stepsPerFrame*SIM.TIME.endStep) = struct('cdata',[],'colormap',[]);

% AXES SETTINGS
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));
ax = gca();
ax.NextPlot = 'replaceChildren';
frame = 1;
for step = 1:SIM.TIME.endStep
    % CHECK IF A FRAME IS TO BE CAPTURED THIS STEP
    if rem(step,stepsPerFrame) ~= 0
        continue
    end
    % STEP DATA
    xData = squeeze(globalTraces(:,1,:,step));
    yData = squeeze(globalTraces(:,2,:,step));
    zData = squeeze(globalTraces(:,3,:,step));
    
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % INITIAL TRACE PLOT
        traceHandle = plot3(ax,xData,yData,zData);
        % UPDATE OBJECT HANDLES
        for entity = 1:SIM.totalObjects
            if isstruct(SIM.OBJECTS(entity).patch)
               % GET THE GLOBAL POSE AT STEP
               [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,entity,step));
               % BUILD MARKER
               markerHandle(entity) = patch(ax,'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + globalMarker(1:3,entity,step)',...
                   'Faces',SIM.OBJECTS(entity).patch.faces,...
                   'FaceColor',SIM.OBJECTS(entity).colour,...
                   'EdgeColor',DATA.figureProperties.EdgeColor,...
                   'EdgeAlpha',DATA.figureProperties.EdgeAlpha,...          
                   'FaceLighting',DATA.figureProperties.FaceLighting,...
                   'LineWidth',DATA.figureProperties.PatchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(entity) = plot3(ax,globalMarker(1,entity,step),globalMarker(2,entity,step),globalMarker(3,entity,step));
                % ASSIGN MARKER PROPERTIES
                markerHandle(entity).Marker = SIM.OBJECTS(entity).symbol;
                markerHandle(entity).MarkerSize = DATA.figureProperties.MarkerSize;
                markerHandle(entity).MarkerFaceColor = SIM.OBJECTS(entity).colour;
                markerHandle(entity).MarkerEdgeColor = DATA.figureProperties.MarkerEdgeColor;
                markerHandle(entity).Color = SIM.OBJECTS(entity).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(entity).LineStyle = DATA.figureProperties.LineStyle;
            traceHandle(entity).LineWidth = DATA.figureProperties.LineWidth;
            traceHandle(entity).Color = SIM.OBJECTS(entity).colour;
            % ANNOTATION WITH THE CURRENT TIME
            clockHandle = annotation('textbox',[0.02 0.05 0.15 0.04],'String','Time:','FitBoxToText','off');
            set(clockHandle,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight);
        end
        % FIGURE PROPERTIES
        title(titleString,'fontweight',DATA.figureProperties.fontWeight,'fontsize',DATA.figureProperties.titleFontSize);
        xlabel('x(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        ylabel('y(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        zlabel('z(m)','fontweight',DATA.figureProperties.fontWeight,'fontSize',DATA.figureProperties.axisFontSize);
        set(ax,'FontSize',DATA.figureProperties.axisFontSize,'fontWeight',DATA.figureProperties.fontWeight,'FontSmoothing','on');
        set(ax,'Color',DATA.figureProperties.axesColor);
%         xlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
%         ylim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
%         zlim([DATA.figureProperties.minPosition,DATA.figureProperties.maxPosition]);
        axis vis3d equal
        set(ax,'outerposition',[0.02 0.1 1 0.88]);
        grid on; box on; hold off;
        legend(markerHandle,DATA.figureProperties.legendEntries,'Location','northeastoutside');
        view([-45 50]);
    else
        % FOR EACH PLOT ENTITY
        for entity = 1:numel(markerHandle)
            if isstruct(SIM.OBJECTS(entity).patch)
                % GET THE GLOBAL POSE AT STEP
                [~,R2] = OMAS_axisTools.quaternionToRotationMatrix(globalMarker(7:end,entity,step));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',SIM.OBJECTS(entity).patch.vertices*R2 + globalMarker(1:3,entity,step)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',globalMarker(1,entity,step));
                set(markerHandle(entity),'YData',globalMarker(2,entity,step));
                set(markerHandle(entity),'ZData',globalMarker(3,entity,step));
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', xData(:,entity));
            set(traceHandle(entity),'YData', yData(:,entity));
            set(traceHandle(entity),'ZData', zData(:,entity));
        end
        % UPDATE TIMESTAMP ANNOTATION
        set(clockHandle,'String',sprintf('Time: %ss',num2str(SIM.TIME.timeVector(1,step))));
    end
    drawnow();
    
    % COLLECT FRAMES
    F(frame) = getframe(figureHandle);
    im = frame2im(F(frame));
    [imind,cm] = rgb2ind(im,256);
    
    % APPEND THE FRAMES TO GIF
    if frame == 1
        imwrite(imind,cm,filePath,'gif','Loopcount',inf,'DelayTime',(1/fps));
    else
        imwrite(imind,cm,filePath,'gif','WriteMode','append','DelayTime',(1/fps));
    end
    frame = frame + 1;                  % Move to next frame
end

% FIGURE COMPLETE
DATA.figureProperties.alignment = DATA.figureProperties.alignment + DATA.figureProperties.spacing;
currentFigure = currentFigure + 1;
end