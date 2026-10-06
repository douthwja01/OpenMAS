% INTERSECT - RANGE
function [intersectFlag] = CheckRangeIntersection(minA,maxA,minB,maxB)
% This function determines the intersection between two 1D
% ranges, and the distance between them.
assert(minA < maxA,'Minimum value must be less than the maximum');
assert(minB < maxB,'Minimum value must be less than the maximum');
% Do either of the ranges on an assumed axis overlap
if maxA <= minB
    intersectFlag = 0;              % They intersect
elseif minA >= maxB
    intersectFlag = 0;              % They intersect
else
    intersectFlag = 1;              % They don't
end

end