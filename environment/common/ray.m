% DEFINE RAY
function [rayObj] = ray(origin,direction,magnitude)
% Define a ray with a given a origin and direction.
assert(numel(origin) == 3 && numel(direction) == 3,'Parameters must be 3D');
if nargin < 3
    magnitude = 1;
end
% NORMALISE THE DIRECTION VECTOR
direction = direction/norm(direction);
% DEFINE THE RAY STRUCTURE
rayObj = struct('origin',origin,...
                'direction',direction,...
                'magnitude',magnitude);
end