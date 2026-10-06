% INTERSECT - SPHERES
function [intersectFlag] = CheckSphereSphereIntersection(positionA,radiusA,positionB,radiusB)
% This uses a simple one dimensional comparison to determine
% the overlap of two spherical volumes.

% COLLISION LOGIC
intersectFlag = 0;
if norm(positionB - positionA) - (radiusA + radiusB) < 0
    intersectFlag = 1;
end
end