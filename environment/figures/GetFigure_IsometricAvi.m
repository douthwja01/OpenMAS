% GET THE 3D TRAJECTORY TRAIS AS A VIDEO [ UPDATED ]
function [currentFigure,figureHandle] = GetFigure_isometricAvi(SIM,objectIndex,DATA,currentFigure)
% This function generates a new 3D trajectory figure be exporting the
% trajectory data to the output directory before generating the figure.

% CONFIGURE THE PLOT ATTRIBUTES
figureHandle = figure('Name','OpenMAS isometric timelapse (AVI)');
titleString = sprintf('Object global trajectories over a period of %ss',num2str(SIM.TIME.endTime));

% MAXIMISE GRAPH SIZE IN WINDOW
setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.08, 0.90, 0.88]);
set(figureHandle,'Position', DATA.figureProperties.windowSettings);        % [x y width height]
set(figureHandle,'Color',DATA.figureProperties.figureColor);               % Background colour 
if DATA.figureProperties.publish
    set(figureHandle,'MenuBar','none');
    set(figureHandle,'ToolBar','none');
    % set(figureHandle,'Visible','off');
end
% CREATE THE AXES IN THE FIGURE
ax = axes(figureHandle);
% ASSUME THE objectIndex is in the same order as the SIM.OBJECTS

% GET THE SIZE OF THE PLOTTED SET
% The data has already been pre-formatted so that the number of states
% aligns with the number of expected frames in the annimation.
% IMPORT THE OBJECT DATA FROM TEMP FILE
load([SIM.outputPath,SIM.systemFile]);
framesToPlot = size(eval(sprintf('objectID%d',SIM.OBJECTS(1).objectID)),2);   % The variable name in the workspace

% PREPARE THE 'AVI' GENERATOR
vidFile = VideoWriter([SIM.outputPath,'isometric']);
vidFile.FrameRate = DATA.figureProperties.fps;                           % Set the videos fps to match the sample rate
vidFile.Quality = 50;
open(vidFile);

% ax.NextPlot = 'replaceChildren';
inclinationAngle = 62;
viewPoint = -10; 
% viewRate = 0.05;
for frame = 1:framesToPlot
    % Must be done on a step by step base to get the annimations correct.
    % BUILD THE FRAME 
    if frame == 1
        hold on;
        % UPDATE OBJECT HANDLES
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT FOR REFERENCE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            objectStates = eval(variableLabel);                              % Evalutate object states        
            % TRACE INITIAL POSITIONS
            traceHandle(entity) = plot3(ax,objectStates(1,frame),objectStates(2,frame),objectStates(3,frame));
            % THE 
            if numel(objectHandle.GEOMETRY.vertices) > 0
               % GET THE GLOBAL POSE AT STEP
               [R_frame] = OMAS_geometry.quaternionToRotationMatrix(objectStates(7:end,frame));
               % BUILD MARKER
               markerHandle(entity) = patch(ax,...
                   'Vertices',objectHandle.GEOMETRY.vertices*R_frame + objectStates(1:3,frame)',...
                   'Faces',objectHandle.GEOMETRY.faces,...
                   'FaceColor',SIM.OBJECTS(entity).colour,...
                   'EdgeColor',DATA.figureProperties.edgeColor,...
                   'EdgeAlpha',DATA.figureProperties.edgeAlpha,...  
                   'FaceLighting',DATA.figureProperties.faceLighting,...
                   'FaceAlpha',DATA.figureProperties.faceAlpha,...
                   'LineWidth',DATA.figureProperties.patchLineWidth);             % Patch properties
            else
                % INITIAL MARKER PLOT
                markerHandle(entity) = plot3(ax,objectStates(1,frame),objectStates(2,frame),objectStates(3,frame));
                % ASSIGN MARKER PROPERTIES
                markerHandle(entity).Marker = SIM.OBJECTS(entity).symbol;
                markerHandle(entity).MarkerSize = DATA.figureProperties.markerSize;
                markerHandle(entity).MarkerFaceColor = SIM.OBJECTS(entity).colour;
                markerHandle(entity).MarkerEdgeColor = DATA.figureProperties.markerEdgeColor;
                markerHandle(entity).Color = SIM.OBJECTS(entity).colour;
            end
            % ASSIGN TRACE PROPERTIES
            traceHandle(entity).LineStyle = DATA.figureProperties.lineStyle;
            traceHandle(entity).LineWidth = DATA.figureProperties.lineWidth;
            traceHandle(entity).Color = SIM.OBJECTS(entity).colour;
        end
        % ANNOTATION WITH THE CURRENT TIME
        clockHandle = annotation('textbox',[0.025 0.025 0.15 0.06],'String','Time:','FitBoxToText','off');
        set(clockHandle,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontWeight',DATA.figureProperties.fontWeight);
        % Title
        title(ax,...
            titleString,...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.titleFontSize,...
            'FontSmoothing','on');
        % X-Label
        xlabel(ax,'x (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Y-Label
        ylabel(ax,'y (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Z-Label
        zlabel(ax,'z (m)',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on');
        % Axes
        set(ax,...
            'TickLabelInterpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName,...
            'Fontweight',DATA.figureProperties.fontWeight,...
            'FontSize',DATA.figureProperties.axisFontSize,...
            'FontSmoothing','on',...
            'Color',DATA.figureProperties.axesColor,...
            'GridLineStyle','--',...
            'GridAlpha',0.25,...
            'GridColor','k');
%         xlim(ax,[DATA.figureProperties.axisMinimums(1),DATA.figureProperties.axisMaximums(1)]);
%         ylim(ax,[DATA.figureProperties.axisMinimums(2),DATA.figureProperties.axisMaximums(2)]);
%         zlim(ax,[DATA.figureProperties.axisMinimums(3),DATA.figureProperties.axisMaximums(3)]);
        axisLimit = (DATA.figureProperties.objectMaximalRadii + DATA.figureProperties.maxAbsPosition);
        xlim(ax,[-axisLimit,axisLimit]);
        ylim(ax,[-axisLimit,axisLimit]);
        zlim(ax,[-axisLimit,axisLimit]);
        axis manual;
        % Legend
        legend(markerHandle,...
            DATA.figureProperties.legendEntries,...
            'Location','northeastoutside',...
            'Interpreter',DATA.figureProperties.interpreter,...
            'fontname',DATA.figureProperties.fontName);
%         set(ax,'OuterPosition', [.1, .2, .9, .6]); % [xLeft, yBottom, width, height]
        view([viewPoint inclinationAngle]);     % Set initial view angle
        grid on; grid minor;
        box on; hold off; 
    else
        % CONTINUED FRAMES
        % FOR EACH PLOT ENTITY
        for entity = 1:SIM.totalObjects
            % GET THE OBJECT FOR REFERENCE
            objectHandle = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entity).objectID};
            % PULL VARIABLE FROM WORKSPACE
            variableLabel = sprintf('objectID%d',SIM.OBJECTS(entity).objectID); % The variable name in the workspace
            objectStates = eval(variableLabel);                                 % Evalutate object states
            % CHECK IF UPDATE IS NECESSARY
            if any(isnan(objectStates(:,frame)))                                % Object is static, freeze its position
                continue
            end
            % HANDLE DIFFERENT REPRESENTATION
            if numel(objectHandle.GEOMETRY.vertices) > 0
                % GET THE GLOBAL POSE AT STEP
                [R_frame] = OMAS_geometry.quaternionToRotationMatrix(objectStates(7:end,frame));
                % BUILD MARKER
                set(markerHandle(entity),'Vertices',objectHandle.GEOMETRY.vertices*R_frame + objectStates(1:3,frame)');
            else
                % UPDATE MARKER DATA
                set(markerHandle(entity),'XData',objectStates(1,frame));
                set(markerHandle(entity),'YData',objectStates(2,frame));
                set(markerHandle(entity),'ZData',objectStates(3,frame));
            end
            % UPDATE TRACE DATA
            if frame <= DATA.figureProperties.tailLength/SIM.TIME.dt
                % TRANSITIONING PERIOD                    
                tailTrace = objectStates(1:3,1:frame);                     % All points upto the step position
            elseif (frame - DATA.figureProperties.tailLength/SIM.TIME.dt) > 0
                % IF THE TRAIL IS NOW A SUBSET
                tailTrace = objectStates(1:3,(frame-(DATA.figureProperties.tailLength/SIM.TIME.dt)):frame);        % All points upto the step position
            else
                tailTrace = objectStates(1:3,frame);        
            end
            % UPDATE TRACE DATA
            set(traceHandle(entity),'XData', tailTrace(1,:));
            set(traceHandle(entity),'YData', tailTrace(2,:));
            set(traceHandle(entity),'ZData', tailTrace(3,:));
        end
    end
    % UPDATE TIMESTAMP ANNOTATION
    set(clockHandle,'String',sprintf('Time: %ss',num2str(frameTimeVector(frame))));
    % FORCE IMAGE WRITE
    drawnow;
    % COLLECT FRAMES    
    writeVideo(vidFile,getframe(figureHandle));       % Write frame to video
    % INCREMENT
%     viewPoint = viewPoint + viewRate;
end
close(vidFile);
currentFigure = currentFigure + 1;
% CLEAN UP
clearvars -except currentFigure figureHandle
end
