% RAY-SPHERE INTERSECTION
function [p_int] = GetRaySphereIntersection(ray,centroid,radius)
% This function calculates the point(s) of intersection between
% a ray and a defined sphere.

% INPUTS:
% ray      - The
% centroid - The center of the defined sphere.
% radius   - The radius of the defined sphere.

% OUTPUT CONTAINER
p_int = [];
tolerance = (1E-3);

% The projection of the separation on the ray direction
t = dot((centroid-ray.origin),ray.direction);

if t < 0
    % Ray is directed away from the sphere
    return
end

%            if t > ray.magnitude
%                % The ray does not extend that far
%                return
%            end

% Rays minimal radius
ySq = t^2 - (norm(centroid-ray.origin))^2;

% Get distance between the centerline and the point on the
% circumference
rSq = radius^2;

% If the radial distance greater than the projection on the radius
if ySq < rSq
    % The magnitude difference
    dt = sqrt(rSq + ySq);
    % First intersection distance
    t1 = t - dt;
    % Require a minimum segment length
    if dt < tolerance
        p_int = (ray.origin + ray.direction*t1)';
        return
    end
    % Second intersection distance
    t2 = t + dt;
    % The second
    p_int = vertcat(p_int,(ray.origin + ray.direction*t2)');
end
end