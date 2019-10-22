%% THE COLLISION DETECTION ROUTINE CHECK (OMAS_collisionDetection.m) %%%%%%
% This function handles the collision detection routines between two
% objects moving through the 3D environment. The collision flag is returned
% based on the hit box assumptions described below.

function [ABcollided] = OMAS_collisionDetection(META_A,geometry_A,META_B,geometry_B,conditionTolerance)

if nargin < 5
    conditionTolerance = 0;
end

% INPUTS:
% META_A     - The global META structure of object A
% geometry_A - The geometry of object A (assumed defined relative to its origin and principle axes)
% META_B     - The global META structure of object B
% geometry_B - The geometry of object B (""."")

% %%%%%%%%%%%%%%%%%% ASSESS THE COLLISION CONDITIONS %%%%%%%%%%%%%%%%%%%%%%
% Objects are given a parameter defining the type of "hitbox" 
% 'none'      - None collidable
% 'spherical' - Spherical assumption
% 'AABB'      - Axis aligned bounding box
% 'OBB'       - Object aligned bounding box

ABcollided = 0;

% NEITHER OBJECTS HAVE A HIT-BOX (Collisions cannot occur)
if META_A.hitBox == OMAS_hitBoxType.none || META_B.hitBox == OMAS_hitBoxType.none
	return 
end   
    
% SIMPLE PROXIMITY CHECK
preliminaryCheck = OMAS_geometry.intersect_spheres(...
    META_A.globalState(1:3),META_A.radius - 0.5*conditionTolerance,...
    META_B.globalState(1:3),META_B.radius - 0.5*conditionTolerance);

if ~preliminaryCheck
    return      % It is not possible for either objects to meet .. return
end

% OBJECT ROTATION MATRICES
R_A = OMAS_geometry.quaternionToRotationMatrix(META_A.globalState(7:10));
R_B = OMAS_geometry.quaternionToRotationMatrix(META_B.globalState(7:10));    

% HIT-BOX CONDITIONS OF OBJECT A
switch META_A.hitBox
    case OMAS_hitBoxType.spherical
        % Behaviour with respect to a sphere representing a object A.
        switch META_B.hitBox
            case OMAS_hitBoxType.spherical
                % Directly intersect its sphere
            	ABcollided = 1;
                return
            case OMAS_hitBoxType.AABB
                % Rotated vertices used defined the box
                rotatedVertices_B = geometry_B.vertices*R_B + META_B.globalState(1:3)';
                % Intersect the AABB: sphere -> cuboid
                ABcollided = OMAS_geometry.intersect_AABB_sphereCuboid(...
                           META_A.globalState(1:3),META_A.radius,...
                           min(rotatedVertices_B),max(rotatedVertices_B));
                return
            case OMAS_hitBoxType.OBB
                % Get the cuboid scaled from the geometry
                cuboidGeometryB = OMAS_graphics.defineCuboid(...
                                min(geometry_B.vertices),max(geometry_B.vertices));
                % Map cuboid to the global space
                rotatedGeometryB = cuboidGeometryB*R_b + META_B.globalState(1:3)';
                % Intersect the OBB: sphere -> cuboid
                ABcollided = OMAS_geometry.intersect_OBB_sphereCuboid(...
                           META_A.globalState(1:3),META_A.radius,...
                           META_B.globalState(1:3),rotatedGeometryB);
                return
            otherwise
                error('Object %d hit box type not recognised.',META_B.objectID); 
        end
    case OMAS_hitBoxType.AABB
        % Behaviour with respect to an AABB representing a object A.
        switch META_B.hitBox
            case OMAS_hitBoxType.spherical
                % Rotated vertices used defined the box
                rotatedVertices_A = geometry_A.vertices*R_A + META_A.globalState(1:3)';
                % Intersect the AABB: sphere -> cuboid
                ABcollided = OMAS_geometry.intersect_AABB_sphereCuboid(...
                           META_B.globalState(1:3),META_B.radius,...
                           min(rotatedVertices_A),max(rotatedVertices_A));
                return
            case OMAS_hitBoxType.AABB
                % Rotated vertices used defined the box
                rotatedVertices_A = geometry_A.vertices*R_A + META_A.globalState(1:3)';
                rotatedVertices_B = geometry_B.vertices*R_B + META_B.globalState(1:3)';
                minA = min(rotatedVertices_A); maxA = max(rotatedVertices_A);
                minB = min(rotatedVertices_B); maxB = max(rotatedVertices_B);
                % Intersect the double AABB: 
                ABcollided = OMAS_geometry.intersect_AABB_cuboids(minA,maxA,minB,maxB);
                return
            case OMAS_hitBoxType.OBB
                % Rotated vertices used defined the box
                rotatedVertices_A = geometry_A.vertices*R_A + META_A.globalState(1:3)';
                % Get the AABB
                cuboidGeometry_A = OMAS_graphics.defineCuboid(min(rotatedVertices_A),...
                                                              max(rotatedVertices_A));
                % Get the OBB
                cuboidGeometry_B = OMAS_graphics.defineCuboid(min(geometry_B.vertices),...
                                                              max(geometry_B.vertices));
                cuboidGeometry_B = cuboidGeometry_B.vertices*R_B + META_B.globalState(1:3)';
                % Intersect the AABB and the OBB
                ABcollided = OMAS_geometry.intersect_OBB_cuboids(...
                           META_A.globalState(1:3),cuboidGeometry_A,...
                           META_B.globalState(1:3),cuboidGeometry_B);
                return
            otherwise
                error('Object %d hit box type not recognised.',META_B.objectID); 
        end
    case OMAS_hitBoxType.OBB
        % Behaviour with respect to an OBB representing a object A.
        switch META_B.hitBox
            case OMAS_hitBoxType.spherical
                % Get the cuboid scaled from the geometry
                cuboidGeometry_A = OMAS_graphics.defineCuboid(...
                                min(geometry_A.vertices),max(geometry_A.vertices));
                % Map cuboid to the global space
                rotatedGeometryA = cuboidGeometry_A*R_A + META_A.globalState(1:3)';
                % Intersect the OBB: sphere -> cuboid
                ABcollided = OMAS_geometry.intersect_OBB_sphereCuboid(...
                           META_B.globalState(1:3),META_B.radius,...
                           META_A.globalState(1:3),rotatedGeometryA);
                return
            case OMAS_hitBoxType.AABB
                % Get the OBB
                cuboidGeometry_A = OMAS_graphics.defineCuboid(min(geometry_A.vertices),...
                                                              max(geometry_A.vertices));
                cuboidGeometry_A = cuboidGeometry_A.vertices*R_A + META_A.globalState(1:3)';
                % Rotated vertices used defined the box
                rotatedVertices_B = geometry_B.vertices*R_B + META_B.globalState(1:3)';
                % Get the AABB
                cuboidGeometry_B = OMAS_graphics.defineCuboid(min(rotatedVertices_B),...
                                                              max(rotatedVertices_B));
                % Intersect the AABB and the OBB
                ABcollided = OMAS_geometry.intersect_OBB_cuboids(...
                           META_B.globalState(1:3),cuboidGeometry_B,...
                           META_A.globalState(1:3),cuboidGeometry_A);
                return
            case OMAS_hitBoxType.OBB
                % Get the OBB A
                cuboidGeometry_A = OMAS_graphics.defineCuboid(min(geometry_A.vertices),...
                                                              max(geometry_A.vertices));
                % Get the OBB B
                cuboidGeometry_B = OMAS_graphics.defineCuboid(min(geometry_B.vertices),...
                                                              max(geometry_B.vertices));
                % Rotated vertices used defined the box
                rotatedVertices_A = cuboidGeometry_A.vertices*R_A + META_A.globalState(1:3)';                                     
                rotatedVertices_B = cuboidGeometry_B.vertices*R_B + META_B.globalState(1:3)';                      
                % Intersect the OBB and the OBB
                ABcollided = OMAS_geometry.intersect_OBB_cuboids(...
                           META_A.globalState(1:3),rotatedVertices_A,...
                           META_B.globalState(1:3),rotatedVertices_B);
                return
            otherwise
                error('Object %d hit box type not recognised.',META_B.objectID); 
        end
    otherwise
        error('Object %d hit box type not recognised.',META_A.objectID);    
end

end