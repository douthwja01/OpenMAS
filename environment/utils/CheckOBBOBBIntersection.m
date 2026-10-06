% INTERSECT - TWO ROTATED CUBOIDS (OBB)
function [intersectFlag,separation_OBB]  = CheckOBBOBBIntersection(centerA,cuboidVerticesA,centerB,cuboidVerticesB)

% RELATIVE POSITION RAYS
ABdirection = (centerB-centerA)/norm(centerB-centerA);
[ABray] = ray(centerA, ABdirection);
[BAray] = ray(centerB,-ABdirection);

% PROJECTIONS OF EACH VERTEX ONTO THE SEPARATION VECTOR
projections_A = zeros(size(cuboidVerticesA,1),1);
for i = 1:size(cuboidVerticesA,1)
    % Its projection on the separation vector
    [projections_A(i)] = GetPointProjectionOnRay(ABray,cuboidVerticesA(i,:)');
end
projections_B = zeros(size(cuboidVerticesB,1),1);
for i = 1:size(cuboidVerticesB,1)
    % Its projection on the separation vector
    [projections_B(i)] = GetPointProjectionOnRay(BAray,cuboidVerticesB(i,:)');
end
% MAXIMAL PROJECTIONS ALONG THE SEPARATION VECTOR
tmax_A = max(projections_A);
tmax_B = max(projections_B);
% CHECK IF THIS EXCEEDS THE SEPARATION OF THE TWO OBJECTS
separation_OBB = norm(centerB-centerA) - (tmax_A + tmax_B);
% EVALUATE INTERSECTION
intersectFlag = 0;
if separation_OBB < 0
    separation_OBB = 0;                                            % Catch negative separations
    intersectFlag = 1;
end
end