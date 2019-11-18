
function [Ry] = R_y(theta)

% Input sanity check
if nargin < 1
    theta = sym('theta_t','real');
end
assert(numel(theta) == 1,'Angle must be a scalar [1x1].');
assert(isa(theta,'sym') || isnumeric(theta),'Angle must be a symbolic or numeric value.');

% Rotate about the y-axis (pitch)
Ry = [cos(theta) 0 sin(theta);
               0 1          0;
     -sin(theta) 0 cos(theta)];
end