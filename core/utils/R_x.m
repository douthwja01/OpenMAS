
function [Rx] = R_x(phi)

% Input sanity check
if nargin < 1
    phi = sym('phi_t','real');
end
assert(numel(phi) == 1,'Angle must be a scalar [1x1].');
assert(isa(phi,'sym') || isnumeric(phi),'Angle must be a symbolic or numeric value.');

% Rotate about the x-axis (roll)
Rx  = [ 1        0         0;
        0 cos(phi) -sin(phi);
        0 sin(phi)  cos(phi)];
end