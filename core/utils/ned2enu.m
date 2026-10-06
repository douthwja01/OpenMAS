% Map vector from ENU to NED
function [v] = ned2enu(v)
% Input sanity check
assert(IsColumn(v,3),'Expecting 3D column vector.');

% Rotate the vector through pi
R_pi = [ 1       0        0;
         0 cos(-pi) -sin(-pi);
         0 sin(-pi)  cos(-pi)];
v = R_pi*v;
end