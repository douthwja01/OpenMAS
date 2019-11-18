
function [Rz] = R_z(psi)

% Input sanity check
if nargin < 1
    psi = sym('psi_t','real');
end
assert(numel(psi) == 1,'Angle must be a scalar [1x1].');
assert(isa(psi,'sym') || isnumeric(psi),'Angle must be a symbolic or numeric value.');

% Rotate about the z-axis (yaw)
Rz = [cos(psi) -sin(psi) 0;
      sin(psi)  cos(psi) 0;
             0         0 1];
end