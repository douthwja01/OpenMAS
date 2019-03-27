% PLANAR PROJECTION OF A VECTOR ON A PLANE
function [v_proj] = GetPlanarProjection(n,v)
% Normalise the normal vector
n = n/norm(n);
% Get the projection on the plane
v_proj = v - (dot(n,v)/norm(n))*n;
end