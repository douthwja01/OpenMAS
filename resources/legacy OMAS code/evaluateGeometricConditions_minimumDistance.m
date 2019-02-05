% GEOMETRIC CONDITIONS BETWEEN TWO OBJECTS

% ASSESS GEOMETRIC CONDITIONS
function [minimumDistance,ABrelativePosition,collisionFlag] = evaluateGeometricConditions(META_A,META_B,objectInd_A,objectInd_B)
% This function calculates the geometric properties necessary for
% evaluating the META event conditions and limits.
% INPUTS:
% META_A - The META.OBJECT for object A
% META_B - The META.OBJECT for object B
% objectInd_A - The class object defining object A
% objectInd_B - The class object defining object B
% OUTPUT:
% minimumDistance  - The scalar distance from surface(point) A to surface(point) B
% relativePosition - The vector between object-central points
% geometricLimit   - The scalar sum of the effective object radii

% Define the new relative global state (x_rel = x_b - x_a)
ABrelativePosition = META_B.globalState(1:3) - META_A.globalState(1:3);    % Relative global position
% Calculate the distance between the two geometries.
relativeABDistance = norm(ABrelativePosition);
unit_relativePosition = ABrelativePosition/relativeABDistance;             % Unit vector from A to B

% DETERMINE BEHAVIOUR BASED ON AVAILABILITY OF GEOMETRY DATA
thisObjHasGeometry   = ~isempty(objectInd_A.GEOMETRY);                     % GEOMETRY property is populated for the object class being updated
secondObjHasGeometry = ~isempty(objectInd_B.GEOMETRY);                     % GEOMETRY property is populated for the object its been evaluated again

CAiterations = 6;
if thisObjHasGeometry && secondObjHasGeometry
    % DESIGN AN INITIAL GUESS TOWARDS THE SECOND OBJECT
    fprintf( 'Both have geometry\n');
    
    fprintf( '%.0f and %.0f\n',objectInd_A.objectID,objectInd_B.objectID);
    
    % GET THE GLOBAL POSE OF THE TWO OBJECTS
    globalPose_A.vertices = objectInd_A.GEOMETRY.vertices*META_A.R + META_A.globalState(1:3)';
    globalPose_A.faces    = objectInd_A.GEOMETRY.faces;
    globalPose_B.vertices = objectInd_B.GEOMETRY.vertices*META_B.R + META_B.globalState(1:3)';
    globalPose_B.faces    = objectInd_B.GEOMETRY.faces;
    
    % INITIAL QUERY POINT
    queryPoint = META_A.globalState(1:3) + unit_relativePosition;
    
    % GET THE MINIMUM SEPARATION FROM A TO B 
    [~,surfacePoint_B] = OMAS_evaluateCollision(globalPose_B,...
        'QueryPoints',queryPoint',...
        'maxDistance',relativeABDistance);                                 % Maximum distance is that to the CG    
     % GET THE MINIMUM SEPARATION FROM B TO A 
    [~,surfacePoint_A] = OMAS_evaluateCollision(globalPose_A,...
        'QueryPoints',surfacePoint_B,...
        'maxDistance',relativeABDistance);                                 % Maximum distance is that to the CG
    
    % ASSESS ONE GEOMETRY AGAINST ANOTHER
    collisionFlag = OMAS_GJK(objectInd_A.GEOMETRY,objectInd_B.GEOMETRY,CAiterations);
    
    % MINIMUM DISTANCE
    minimumDistance = norm(surfacePoint_B - surfacePoint_A);
    if collisionFlag
        minimumDistance = 0;
    end
elseif  thisObjHasGeometry && ~secondObjHasGeometry
    % ASSESS ONE GEOMETRY AGAINST A SPHERICAL CONSTRAINT
    fprintf( 'A does, B does not \n');
    
    % GET THE GLOBAL POSE OF THE TWO OBJECTS
    globalPose_A.vertices = objectInd_A.GEOMETRY.vertices*META_A.R + META_A.globalState(1:3)';
    globalPose_A.faces    = objectInd_A.GEOMETRY.faces;
    
    % INITIAL QUERY POINT
    radialPoint_B = META_B.globalState(1:3) - META_B.radius*unit_relativePosition;
    
    % GET THE MINIMUM SEPARATION FROM A TO B 
    [~,surfacePoint_A] = OMAS_evaluateCollision(globalPose_A,...
        'QueryPoints',radialPoint_B',...
        'maxDistance',relativeABDistance);                                 % Maximum distance is that to the CG
    
    % SPHERICAL APPROXIMATION OF BODY B
    [X,Y,Z] = sphere(20);
    X = X.*META_B.radius + META_B.globalState(1);
    Y = Y.*META_B.radius + META_B.globalState(2);
    Z = Z.*META_B.radius + META_B.globalState(3);
    geometryAssumption_B = surf2patch(X,Y,Z);
    
    % ASSESS ONE GEOMETRY AGAINST ANOTHER
    collisionFlag = OMAS_GJK(objectInd_A.GEOMETRY,geometryAssumption_B,CAiterations);
    
    % MINIMUM DISTANCE
    minimumDistance = norm(radialPoint_B - surfacePoint_A);
    if collisionFlag
        minimumDistance = 0;
    end
elseif ~thisObjHasGeometry && secondObjHasGeometry
    % ASSESS ONE GEOMETRY AGAINST A SPHERICAL CONSTRAINT
    fprintf( 'A does not, B Does. \n');
    
    % GET THE GLOBAL POSE OF THE TWO OBJECTS
    globalPose_B.vertices = objectInd_B.GEOMETRY.vertices*META_B.R + META_B.globalState(1:3)';
    globalPose_B.faces    = objectInd_B.GEOMETRY.faces;
    
    % INITIAL QUERY POINT
    radialPoint_A = META_A.globalState(1:3) + META_A.radius*unit_relativePosition;
    
    % GET THE MINIMUM SEPARATION FROM A TO B 
    [~,surfacePoint_B] = OMAS_evaluateCollision(globalPose_B,...
        'QueryPoints',radialPoint_A',...
        'maxDistance',relativeABDistance);                                 % Maximum distance is that to the CG
    
    % SPHERICAL APPROXIMATION OF BODY A
    [X,Y,Z] = sphere(20);
    X = X.*META_A.radius + META_A.globalState(1);
    Y = Y.*META_A.radius + META_A.globalState(2);
    Z = Z.*META_A.radius + META_A.globalState(3);
    geometryAssumption_A = surf2patch(X,Y,Z);
    
    % ASSESS ONE GEOMETRY AGAINST ANOTHER
    collisionFlag = OMAS_GJK(objectInd_A.GEOMETRY,geometryAssumption_A,CAiterations);
    
    % MINIMUM DISTANCE
    minimumDistance = norm(surfacePoint_B - radialPoint_A);
    if collisionFlag
        minimumDistance = 0;
    end
    
else  
    % ASSESS ONE SPHERICAL CONSTRAINT AGAINST ANOTHER SPHERICAL CONSTRAINT
    fprintf( 'Neither A or B have geometries \n');
    
    minimumDistance = relativeABDistance - (META_A.radius + META_B.radius);
    collisionFlag = 0;
    if minimumDistance < 0
        collisionFlag = 1;
    end
end


% THE MINIMUM ALLOWABLE DISTANCE
% geometricLimit  = norm(relativePosition) - minimumDistance;

% Enabled = 0;
% if META_A.objectID == 1 && Enabled
%     % FIGURE PREPARATION
%     ax = gca;
%     hold on;
%     grid on;
%     axis square equal;
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     zlabel('Z (m)');
%     view([-120 35]);
% 
%     % PLOT THE FIRST AGENT VERTEX DATA
%     if ~isempty(META_A.geometry)
%         patch(ax,'Vertices',META_A.geometry.vertices*META_A.R + META_A.globalState(1:3)',...
%             'Faces',META_A.geometry.faces,...
%             'FaceColor','r',...
%             'EdgeColor','k',...
%             'EdgeAlpha',0.2,...
%             'FaceLighting','gouraud',...
%             'LineWidth',1);
%     end
%     [sphereHandle] = OMAS_axisTools.drawSphere(gcf,META_A.globalState(1:3),META_A.radius);
% 
%     % PLOT THE SECOND AGENT VERTEX DATA
%     if ~isempty(META_B.geometry)
%         patch(ax,'Vertices',META_B.geometry.vertices*META_B.R + META_B.globalState(1:3)',...
%             'Faces',META_B.geometry.faces,...
%             'FaceColor','g',...
%             'EdgeColor','k',...
%             'EdgeAlpha',0.2,...
%             'FaceLighting','gouraud',...
%             'LineWidth',1);
%     end
%     [sphereHandle] = OMAS_axisTools.drawSphere(gcf,META_B.globalState(1:3),META_B.radius);
% 
%     % PLOT THE POINT EVALUATED
%     scatter3(ax,pointB(1),pointB(2),pointB(3),'k*'); % Evaluated point on the geometries
%     scatter3(ax,pointA(1),pointA(2),pointA(3),'m');
% 
%     projection = pointB - pointA;
%     q = quiver3(pointA(1),pointA(2),pointA(3),...
%         projection(1),projection(2),projection(3),'r');
%     q.AutoScaleFactor = 1;
% end
end