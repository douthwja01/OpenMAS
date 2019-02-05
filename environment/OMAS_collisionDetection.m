%% THE COLLISION DETECTION ROUTINE CHECK (OMAS_collisionDetection.m) %%%%%%
% This function handles the collision detection routines between two
% objects moving through the 3D environment. The collision flag is returned
% based on the hit box assumptions described below.

function [ABCollided] = OMAS_collisionDetection(META_A,geometry_A,META_B,geometry_B)

% INPUTS:
% META_A     - The global META structure of object A
% geometry_A - The geometry of object A (assumed defined relative to its origin and principle axes)
% META_B     - The global META structure of object B
% geometry_B - The geometry of object B (""."")

% COLLISION DETECTION (BASED ON OBJECT ASSUMPTIONS)
% The separation check is what prevents us from using early returns.
% We have the following scenarios:
% - Object A is a polygon ; Object B is a polygon
% - Object A is a polygon ; Object B is a point
% - Object A is a point   ; Object B is a polygon
% - Object A is a point   ; Object B is a point

% DETERMINE BEHAVIOUR BASED ON AVAILABILITY OF GEOMETRY DATA
objectA_isPoint = size(geometry_A.vertices,1) < 1;                         % GEOMETRY property is populated for the object class being updated
objectB_isPoint = size(geometry_B.vertices,1) < 1;                         % GEOMETRY property is populated for the object its been evaluated again

% ///////////// ASSESS THE COLLISION SEPARATION & COLLISION /////////////// 
if ~objectA_isPoint && ~objectB_isPoint
    % COMPARING COMPLEX OBSTACLE TO COMPLEX OBSTACLE
    % CONSTRAINING VOLUME A
    dimMin = min(geometry_A.vertices,[],1); 
    dimMax = max(geometry_A.vertices,[],1);
    [cuboidA] = OMAS_graphics.defineCuboid(dimMin,dimMax);         % A's cuboid volume
    % CONSTRAINING VOLUME B
    dimMin = min(geometry_B.vertices,[],1); 
    dimMax = max(geometry_B.vertices,[],1);
    [cuboidB] = OMAS_graphics.defineCuboid(dimMin,dimMax);         % B's cuboid volume
    % MOVE AND ORIENTATE THE CUBOID
    ABrelativePosition = META_B.globalState(1:3,1) - META_A.globalState(1:3,1);
    relativeVertices = cuboidB.vertices*(META_A.R*META_B.R) + ABrelativePosition';
    % COMPARE TWO ORIENTATED CUBOIDS 
    [ABCollided,~] = OMAS_geometry.intersect_OBB_cuboids(zeros(3,1),...
                                                         cuboidA.vertices,...
                                                         ABrelativePosition,...
                                                         relativeVertices);
elseif ~objectA_isPoint
    % COMPARING COMPLEX OBSTACLE TO (AGENT,WAYPOINT,SIMPLE OBSTACLES)
    % CONSTRAINING VOLUME A
    dimMinA = min(geometry_A.vertices,[],1); 
    dimMaxA = max(geometry_A.vertices,[],1);
    [cuboidA] = OMAS_graphics.defineCuboid(dimMinA,dimMaxA);
    % MOVE AND ORIENTATE THE CUBOID
    cuboidVerticesA = cuboidA.vertices*META_A.R + META_A.globalState(1:3)'; % Map a cuboid to the global space of A
    % COMPARE ORIENTATED CUBOID WITH SPHERE
    [ABCollided,~] = OMAS_geometry.intersect_OBB_sphereCuboid(META_B.globalState(1:3),...
                                                              META_B.radius,...
                                                              META_A.globalState(1:3),...
                                                              cuboidVerticesA);
elseif ~objectB_isPoint
    % COMPARING COMPLEX OBSTACLE TO (AGENT,WAYPOINT,SIMPLE OBSTACLES)
    dimMinB = min(geometry_B.vertices,[],1); 
    dimMaxB = max(geometry_B.vertices,[],1);
    [cuboidB] = OMAS_graphics.defineCuboid(dimMinB,dimMaxB);
    cuboidVerticesB = cuboidB.vertices*META_B.R + META_B.globalState(1:3)'; % Map A into the frame of B
    % COMPARE ORIENTATED CUBOID WITH SPHERE
    [ABCollided,~] = OMAS_geometry.intersect_OBB_sphereCuboid(META_A.globalState(1:3),...
                                                              META_A.radius,...
                                                              META_B.globalState(1:3),...
                                                              cuboidVerticesB);
else
    % NO COMPLEX COMPARISON REQUIRED: ASSUME SPHERES
    [ABCollided,~] = OMAS_geometry.intersect_spheres(...
        META_A.globalState(1:3),META_A.radius,...
        META_B.globalState(1:3),META_B.radius);
    return  % Separation is trivial
end

end


