% SIGNED ANGLE BETWEEN TO VECTORS
function [planarAngle,isRightLogical] = GetSignedPlanarAngle(n,u,v)
% This function computes the signed angle between two vectors
% with respect to a provided normal.

% RENORMALISE THE INPUTS
n = unit(n);
u = unit(u);
v = unit(v);
% DETERMINE THE DIRECTION (CW+ve with respect to the normal)
isRightLogical = -sign(dot(cross(u,v),n));
% DETERMINE THE MAGNITUDE
[planarAngle] = GetPlanarAngle(u,v);
% COMBINE
planarAngle = isRightLogical*planarAngle;
end