%% EVALUATE THE RELATIONSHIPS BETWEEN TWO 3D GEOMETRIES (OMAS_geometricRelations.m)
% This function is intended to generate the properties defining the
% geometric relationship between two objects in 3D space.

function [ABrelativePosition,ABSeparation,ABCollided] = OMAS_geometricRelations(META_A,META_B,object_A,object_B)
% INPUTS:
% META_A   - The META structure for object_A
% META_B   - The META strucutre for object_B
% object_A - The class instance of object_A 
% object_B - The class instance of object_B

% GET THE RELATIVE POSITION OF THE TWO OBJECTS (CG-CG)
ABrelativePosition = META_B.globalState(1:3,1) - META_A.globalState(1:3,1);

%% /////////////////// COLLISION DETECTION ROUTINE ////////////////////////
% USE THE AABB APPROXIMATION TO DETERMINE COLLISION
[ ABCollided ] = collisionCheck_AABB(META_A,META_B,object_A,object_B);     % Check for AABB violation
fprintf('\tCheck A\t%s & %s collision: %.0f\n',META_A.name,META_B.name,ABCollided);
if ABCollided
    % USE THE OBB  APPROXIMATION TO DETERMINE COLLISION
    [ ABCollided ] = collisionCheck_OBB(META_A,META_B,object_A,object_B);  % If possible collision, check OBB violation
%     fprintf('\tCheck B\t%s & %s collision: %.0f\n',META_A.name,META_B.name,ABCollided);
end

% DETERMINE BEHAVIOUR BASED ON AVAILABILITY OF GEOMETRY DATA
objectA_hasGeometry = ~isempty(object_A.GEOMETRY);                         % GEOMETRY property is populated for the object class being updated
objectB_hasGeometry = ~isempty(object_B.GEOMETRY);                         % GEOMETRY property is populated for the object its been evaluated again

%% "SIMPLE OBJECT" APPROXIMATION (NEITHER OBJECT HAS GEOMETRY)
if ~objectA_hasGeometry && ~objectB_hasGeometry
    % THE SEPARATION OF TWO SPHERE APPROXIMATIONS
    ABSeparation = norm(ABrelativePosition) - (META_A.radius + META_B.radius);
    return
end

%% DECLARE COMPLEX GEOMETRY ASSUMPTIONS
if objectA_hasGeometry && objectB_hasGeometry
    % DEFINE THE POSE OF A
    globalPose_A.vertices = object_A.GEOMETRY.vertices*META_A.R + META_A.globalState(1:3,1)';
    globalPose_A.faces    = object_A.GEOMETRY.faces;
    % DEFINE THE POSE OF B
    globalPose_B.vertices = object_B.GEOMETRY.vertices*META_B.R + META_B.globalState(1:3,1)';
    globalPose_B.faces    = object_B.GEOMETRY.faces;
elseif objectA_hasGeometry && ~objectB_hasGeometry
    % DEFINE THE POSE OF A
    globalPose_A.vertices = object_A.GEOMETRY.vertices*META_A.R + META_A.globalState(1:3,1)';
    globalPose_A.faces    = object_A.GEOMETRY.faces;
    % SPHERICAL APPROXIMATION OF BODY B
    [X,Y,Z] = sphere(10);
    X = X.*META_B.radius + META_B.globalState(1);
    Y = Y.*META_B.radius + META_B.globalState(2);
    Z = Z.*META_B.radius + META_B.globalState(3);
    globalPose_B = surf2patch(X,Y,Z,'triangles');                                      % Approximate B as a sphere
else
    % DEFINE THE POSE OF B
    globalPose_B.vertices = object_B.GEOMETRY.vertices*META_B.R + META_B.globalState(1:3,1)';
    globalPose_B.faces    = object_B.GEOMETRY.faces;
    % SPHERICAL APPROXIMATION OF BODY A
    [X,Y,Z] = sphere(10);
    X = X.*META_A.radius + META_A.globalState(1);
    Y = Y.*META_A.radius + META_A.globalState(2);
    Z = Z.*META_A.radius + META_A.globalState(3);
    globalPose_A = surf2patch(X,Y,Z,'triangles');                                      % Approximate A as a sphere
end

%% WE NOW COMPARE THE TWO GEOMETRIES
% We now need to evaluate the minimum distance between the geometries, in
% addition to the collision condition.
queryPoint0 = META_A.globalState(1:3,1) + META_A.radius*ABrelativePosition/norm(ABrelativePosition);
searchSpaceLimit = norm(ABrelativePosition);
% GET THE MINIMUM SEPARATION FROM A TO B
[~,surfacePoint_B] = OMAS_evaluateCollision(globalPose_B,...
    'QueryPoints',queryPoint0',...
    'maxDistance',searchSpaceLimit);                                       % Maximum distance is that to the CG
% GET THE MINIMUM SEPARATION FROM B TO A
[~,surfacePoint_A] = OMAS_evaluateCollision(globalPose_A,...
    'QueryPoints',surfacePoint_B,...
    'maxDistance',searchSpaceLimit);                                       % Maximum distance is that to the CG

% THE SURFACE TO SURFACE VECTOR
minimumVector = surfacePoint_B - surfacePoint_A;
ABSeparation = norm(minimumVector);
if ABCollided 
    ABSeparation = 0;
end
end

% COLLISION DETECTION: ORIENTATED BOUNDING BOX (OBB) 
function [collisionFlag,separation] = collisionCheck_OBB(META_A,META_B,object_A,object_B)
% This function computes the collision detection for two geometries by
% approximating them as two object-aligned cuboids

% DETERMINE BEHAVIOUR BASED ON AVAILABILITY OF GEOMETRY DATA
objectA_hasGeometry = ~isempty(object_A.GEOMETRY);                         % GEOMETRY property is populated for the object class being updated
objectB_hasGeometry = ~isempty(object_B.GEOMETRY);                         % GEOMETRY property is populated for the object its been evaluated again

% COMPUTE THE X,Y,Z INTERVALS IN THE OBJECT SPACE OF A
if ~objectA_hasGeometry
   % MIN/MAX COORDINATES FOR THE OBJECT
    dimMinA = -ones(3,1)*META_A.radius/1.7321;
    dimMaxA =  ones(3,1)*META_A.radius/1.7321;                             % Create a cuboid within the radius
else
    % GET THE VERTEX MAXIMALS/MINIMALS
    vertices = object_A.GEOMETRY.vertices;
    dimMinA = min(vertices,[],1);
    dimMaxA = max(vertices,[],1);
end

if ~objectB_hasGeometry
    % REPRESENT THE OBJECT IN THE BODY AXES OF A
    dimMinB = -ones(3,1)*META_B.radius/1.7321;
    dimMaxB =  ones(3,1)*META_B.radius/1.7321;                             % Create a cuboid within the radius
    
else
    % ROTATE THE GEOMETRY OF B INTO A 
    vertices = object_B.GEOMETRY.vertices;
    dimMinB = min(vertices,[],1);
    dimMaxB = max(vertices,[],1);
end

% DECLARE PARAMETERS
centerA = zeros(3,1);
centerB = META_B.globalState(1:3,1) - META_A.globalState(1:3,1);
norm_AB = norm(centerB);
unit_AB = centerB/norm_AB; 
% THE ROTATION OF OBJECT B WITH RESPECT TO A'S FRAME
R_B = (META_A.R*META_B.R);

% DEFINE THE RAY BETWEEN CENTROIDS
ABray = OMAS_geometry.definedRay(centerA,unit_AB);
BAray = OMAS_geometry.definedRay(centerB,-unit_AB);
% DEFINE VERTEX DATA
cuboidVerticesA = OMAS_geometry.defineCuboid(dimMinA,dimMaxA);
cuboidVerticesB = OMAS_geometry.defineCuboid(dimMinB,dimMaxB);
cuboidVerticesB = cuboidVerticesB*R_B + centerB';               % Rotate and translate the vertices

% PROJECTIONS OF EACH VERTEX ONTO THE SEPARATION VECTOR
projections_A = zeros(size(cuboidVerticesA,1),1);
for i = 1:size(cuboidVerticesA,1)
    queryPoint = cuboidVerticesA(i,:)' - centerA;
    % Its projection on the separation vector
    [projections_A(i)] = OMAS_geometry.pointProjectionOnRay(ABray,queryPoint);
end
projections_B = zeros(size(cuboidVerticesB,1),1);
for i = 1:size(cuboidVerticesB,1)
    queryPoint = cuboidVerticesB(i,:)';
    % Its projection on the separation vector
    [projections_B(i)] = OMAS_geometry.pointProjectionOnRay(BAray,queryPoint);
end
% MAXIMAL PROJECTIONS ALONG THE SEPARATION VECTOR
tmax_A = max(projections_A);
tmax_B = max(projections_B);
% CHECK IF THIS EXCEEDS THE SEPARATION OF THE TWO OBJECTS
separation = norm_AB - (tmax_A + tmax_B);

% DEBUG: PLOT THE OBB SCENE
Enabled = 0;
if META_A.objectID == 1 && Enabled
    figure(1)
    axis vis3d; hold on; grid on;
    scatter3(cuboidVerticesA(:,1),cuboidVerticesA(:,2),cuboidVerticesA(:,3),'ro');
    scatter3(cuboidVerticesB(:,1),cuboidVerticesB(:,2),cuboidVerticesB(:,3),'bo');
    plot3([centerA(1);centerB(1)],[centerA(2);centerB(2)],[centerA(3);centerB(3)],'k');
    
    % PLOT THE MAXIMUM PROJECTION TOWARD EACH OTHER
    q = quiver3(ABray.origin(1),ABray.origin(2),ABray.origin(3),...
        tmax_A*ABray.direction(1),tmax_A*ABray.direction(2),tmax_A*ABray.direction(3),'g');
    q.AutoScaleFactor = 1;
    
    q = quiver3(BAray.origin(1),BAray.origin(2),BAray.origin(3),...
        tmax_B*BAray.direction(1),tmax_B*BAray.direction(2),tmax_B*BAray.direction(3),'g');
    q.AutoScaleFactor = 1;
    cla;
end

% DETETMINE IF THERES BEEN A COLLISION
collisionFlag = 0;
if separation < 0
    separation = 0;
    collisionFlag = 1;
    return
end
end
% COLLISION DETECTION: AXIS-ALIGNED BOUNDING BOX (AABB)
function [collisionFlag] = collisionCheck_AABB(META_A,META_B,object_A,object_B)

% DETERMINE BEHAVIOUR BASED ON AVAILABILITY OF GEOMETRY DATA
objectA_hasGeometry = ~isempty(object_A.GEOMETRY);                         % GEOMETRY property is populated for the object class being updated
objectB_hasGeometry = ~isempty(object_B.GEOMETRY);                         % GEOMETRY property is populated for the object its been evaluated again

% COMPUTE THE X,Y,Z INTERVALS EITHER OBJECT BELONG TO
if ~objectA_hasGeometry
   % MIN/MAX COORDINATES FOR THE OBJECT
    dimMinA = META_A.globalState(1:3) - META_A.radius/1.7321;
    dimMaxA = META_A.globalState(1:3) + META_A.radius/1.7321;              % Dimensions of a cube internal to the sphere
else
    vertices = object_A.GEOMETRY.vertices*META_A.R + META_A.globalState(1:3,1)'; % 
    dimMinA = min(vertices,[],1);
    dimMaxA = max(vertices,[],1);
end

if ~objectB_hasGeometry
    % MIN/MAX COORDINATES FOR THE OBJECT
    dimMinB = META_B.globalState(1:3) - META_B.radius/1.7321;
    dimMaxB = META_B.globalState(1:3) + META_B.radius/1.7321;
else
    vertices = object_B.GEOMETRY.vertices*META_B.R + META_B.globalState(1:3,1)';
    dimMinB = min(vertices,[],1);
    dimMaxB = max(vertices,[],1);
end
% INTERSECT THE TWO CUBOIDS
[collisionFlag] = OMAS_geometry.intersect_cuboids(dimMinA,dimMaxA,dimMinB,dimMaxB);

end
% COLLISION DETECTION: SPHERICAL APPROXIMATION
function [collisionFlag] = collisionCheck_radial(META_A,META_B,object_A,object_B)
% GLOBAL SEPARATION
distance = norm(META_B.globalState(1:3,1) - META_A.globalState(1:3,1));
% COMBINED GEOMETRIC RADII
combinedRadius = META_A.radius + META_B.radius;                            % The two radii already include all geometric points
% COLLISION LOGIC
collisionFlag = 0;
if distance < combinedRadius
    collisionFlag = 1;
end
end
