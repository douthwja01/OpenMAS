% INTERSECT - SPHERE AND ROTATED CUBOID (OBB)
function [intersectFlag,separation_OBB]  = CheckOBBSphereIntersection(centerA,radiusA,centerB,cuboidVerticesB)
% This function computes the intersection between a sphere and a
% rotated bounding box (OBB). This function assumes that the
% problem is being resolved in the axes of A (i.e) 'R' is the
% rotation matrix taking the cube from its own frame to the
% axes of A.

% RELATIVE POSITION RAY
[BAray] = ray(centerB,(centerA-centerB));       % Vector between
% If we define A to be sphere, its projection on the separation
% vector will always be its radius
projection_A = radiusA;
% GET THE PROJECTION OF ALL THE VERTICES ON THE SEPARATING VECTOR
projection_B = zeros(size(cuboidVerticesB,1),1);
for i = 1:size(cuboidVerticesB,1)
    % Its projection on the separation vector
    [projection_B(i)] = GetPointProjectionOnRay(BAray,cuboidVerticesB(i,:)');
end
% THE MAXIMAL PROJECTION TOWARDS THE SPHERE
projection_B = max(projection_B);
% CHECK IF THIS EXCEEDS THE SEPARATION OF THE TWO OBJECTS
separation_OBB = norm(centerB-centerA) - (projection_A + projection_B);
% EVALUATE INTERSECTION
intersectFlag = 0;
if separation_OBB < 0
    separation_OBB = 0;
    intersectFlag = 1;
end
end