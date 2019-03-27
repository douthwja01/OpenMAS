% GET THE INTERSECTION BETWEEN TWO LINES IN 3D
function [pa,pb,isSuccessful] = GetRayRayIntersection3D(p1,dp1,p2,dp2)
%     /*
%        Calculate the line segment PaPb that is the shortest route between
%        two lines P1P2 and P3P4. Calculate also the values of mua and mub where
%           Pa = P1 + mua (P2 - P1)
%           Pb = P3 + mub (P4 - P3)
%        Return FALSE if no solution exists.
%     */

% SOME SUFFICIENTLY SMALL VALUE FOR INTERSECTION
EPS = 1E-9;
isSuccessful = logical(false);
pa = NaN(3,1);
pb = NaN(3,1);

% GET THE DIRECTIONS
p13 = zeros(3,1);
p13(1) = p1(1) - p2(1);
p13(2) = p1(2) - p2(2);
p13(3) = p1(3) - p2(3);

%             if (abs(dP3(1)) < EPS && abs(dP3(2)) < EPS && abs(dP3(3)) < EPS)
%                 return
%             end
%
%             if (abs(dP1(1)) < EPS && abs(dP1(2)) < EPS && abs(dP1(3)) < EPS)
%                 return
%             end

d1343 = p13(1) * dp2(1) + p13(2) * dp2(2) + p13(3) * dp2(3);
d4321 = dp2(1) * dp1(1) + dp2(2) * dp1(2) + dp2(3) * dp1(3);
d1321 = p13(1) * dp1(1) + p13(2) * dp1(2) + p13(3) * dp1(3);
d4343 = dp2(1) * dp2(1) + dp2(2) * dp2(2) + dp2(3) * dp2(3);
d2121 = dp1(1) * dp1(1) + dp1(2) * dp1(2) + dp1(3) * dp1(3);

denom = d2121 * d4343 - d4321 * d4321;
if (abs(denom) < EPS)
    return
end

numer = d1343 * d4321 - d1321 * d4343;

mua = numer / denom;
mub = (d1343 + d4321 * mua) / d4343;

% BUILD INTERSECTION COORDINATES
pa(1) = p1(1) + mua * dp1(1);
pa(2) = p1(2) + mua * dp1(2);
pa(3) = p1(3) + mua * dp1(3);
pb(1) = p2(1) + mub * dp2(1);
pb(2) = p2(2) + mub * dp2(2);
pb(3) = p2(3) + mub * dp2(3);
% IS SUCCESSFUL
isSuccessful = logical(true);
end