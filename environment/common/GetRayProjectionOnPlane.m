% PROJECTION OF A RAY ON A PLANE DEFINED BY A NORMAL
function [v_prj] = GetRayProjectionOnPlane(ray,planeNormal)
% This function computes the projection of a vector on a plane
% defined by its normal vector.

% Ensure normal vector is the normal
planeNormal = planeNormal/norm(planeNormal);
% Scale the ray to its maximum length
v_prj = ray.magnitude*ray.direction - dot(ray.direction,planeNormal)*planeNormal;

end