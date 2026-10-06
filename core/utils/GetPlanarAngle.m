% PLANAR ANGLE BETWEEN TWO VECTORS
function [planarAngle] = GetPlanarAngle(u,v)
planarAngle = acos(dot(u,v)/(norm(u)*norm(v)));
end