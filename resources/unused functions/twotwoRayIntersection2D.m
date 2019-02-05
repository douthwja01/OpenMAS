% TWO LINE INTERSECTIONS (FROM agent_VO)
function [p_inter,isSuccessful] = twotwoRayIntersection2D(P1,dP1,P2,dP2)
% Find the intersection point between two 2D vectors. This
% function isnt' interested if vertices are infront or behind
% of the starting point.
% INPUTS:
% - P1,P2   - The ray defining points.
% - dP1,dP2 - The ray unit directions.
% OUTPUTS:
% - p_inter - The 2D intersection point.

assert(numel(P1) == 2,'Input must be 2D');
assert(numel(P2) == 2,'Input must be 2D');

% SOME SUFFICIENTLY SMALL VALUE FOR INTERSECTION
isSuccessful = logical(false);   % Default to no intersection
p_inter = NaN(2,1);              % Default to no intersection

% THE 2D DETERMININANT
div = dP1(2)*dP2(1) - dP1(1)*dP2(2);
if div == 0
    disp('Lines are parallel');
    return % Lines are parallel
end

% THE SCALAR PROJECTIONS
mua = (dP2(1)*(P2(2) - P1(2)) + dP2(2)*(P1(1) - P2(1))) / div;
mub = (dP1(1)*(P2(2) - P1(2)) + dP1(2)*(P1(1) - P2(1))) / div;

% POINTS MUST BE THE RESULT OF A POSITIVE INCREMENT OF THE VECTOR GRADIENT
% (i.e, in the correct direction)
if mua < 0 || mub < 0   % Intersections only occur in the direction of the vector
    return              % Lines do not intersect
end

% THE INTERSECTION POINT
p_inter = P1 + mua*dP1;
isSuccessful = logical(true);
end