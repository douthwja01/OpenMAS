% PLOT SCENARIO FROM INITIALISED OBJECT-INDEX
function [figureHandle] = plotObjectIndex(obj,objectIndex)
% This function is designed to plot the configured objectIndex
% using their global properties and thier simulation
% object-types

global plotnum

% DETERMINE PLOT PROPERTIES
if ~exist('plotnum','var') || isempty(plotnum)
    plotnum = 1;    % Default to first plot
end

% GENERATE THE FIGURE
figureHandle = figure(plotnum);
axis('equal');
xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
hold on; grid on;
set(gca,'FontSize',12,'fontWeight','bold');

% OBJECT COUNTERS
agents = 0; obstacles = 0; waypoints = 0;
for index = 1:numel(objectIndex)
    % GET THE OBSTACLES GLOBAL PARAMETERS
    objectName       = objectIndex{index}.name;
    objectRadius     = objectIndex{index}.VIRTUAL.radius;
    objectType       = objectIndex{index}.VIRTUAL.type;
    globalPosition   = objectIndex{index}.VIRTUAL.globalPosition;
    globalVelocity   = objectIndex{index}.VIRTUAL.globalVelocity;
    globalQuaternion = objectIndex{index}.VIRTUAL.quaternion;
    % GET THE ROTATION MATRIX GOING FROM BODY-GLOBAL
    [R_q] = OMAS_axisTools.quaternionToRotationMatrix(globalQuaternion);
    % PLOT THE OBJECT ORIENTATION TRIAD
    colours = {'r','g','b'};
    % REPLOT TRIAD
    for j = 1:size(obj.globalTriad,2)
        % GET THE ROTATED TRIAD AXES
        rotatedVector = R_q'*obj.globalTriad(:,j);
        % DRAW THE LOCAL TRIAD
        quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
            rotatedVector(1),rotatedVector(2),rotatedVector(3),colours{j});
        quiv.AutoScaleFactor = 1;
    end
    % PLOT THE POSITIONS
    scatter3(globalPosition(1),globalPosition(2),globalPosition(3),'r');
    % PLOT THE VELOCITY VECTORS
    quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
        globalVelocity(1),globalVelocity(2),globalVelocity(3),'c');
    quiv.AutoScaleFactor = 1;
    % ADD REPRESENTATIVE SPHERES
    switch objectType
        case OMAS_objectType.agent
            % IS AGENT
            [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
            objectColour = 'b';
            agents = agents + 1;
        case OMAS_objectType.obstacle
            % IS OBSTACLE
            [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
            objectColour = 'r';
            obstacles = obstacles + 1;
        case OMAS_objectType.waypoint
            % IS WAYPOINT
            [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
            objectColour = 'g';
            waypoints = waypoints + 1;
        otherwise
            error('[SCENARIO]\tObject type not recognised');
    end
    sphZone = mesh(X,Y,Z);
    set(sphZone,'facealpha',0.25,...
        'FaceColor',objectColour,...
        'LineWidth',0.1,...
        'EdgeAlpha',0.1);
    % ADD ANNOTATION
    annotationText = sprintf(' \t%s [ID:%s]',objectName,num2str(objectIndex{index}.objectID));
    text(globalPosition(1),globalPosition(2),globalPosition(3),char(annotationText));
end
% ADD TITLE
titleStr = sprintf('Test scenario: %s agents, %s obstacles and %s waypoints.',num2str(agents),num2str(obstacles),num2str(waypoints));
title(titleStr);
hold off;
end