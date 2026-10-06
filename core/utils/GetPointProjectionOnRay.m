% PROJECTION OF A POINT ON A VECTOR
function [v_mag] = GetPointProjectionOnRay(ray,point)
% Calcuate the projection magnitude(distance)
v_mag = dot(ray.direction,(point - ray.origin));
end