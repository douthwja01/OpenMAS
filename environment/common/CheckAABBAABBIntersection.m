% INTERSECT - AXIS ALIGNED CUBOIDS (AABB)
function [intersectFlag] = CheckAABBAABBIntersection(minA,maxA,minB,maxB)
% This function determines if two 3D cuboids currently intersect
intersectFlag = 1;
iMatrix = zeros(numel(minA),2);
for dim = 1:numel(minA)
    % INTERSECTING TWO RANGES, EACH DIMENSION AT A TIME
    iMatrix(dim,1) = CheckRangeIntersection(minA(dim),maxA(dim),minB(dim),maxB(dim));
    if ~iMatrix(dim,1)
        intersectFlag = 0;
    end
end
end