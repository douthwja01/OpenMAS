% Map vector from ENU to NED
function [v] = enu2ned(v)
% Input sanity check
assert(IsColumn(v,3),'Expecting 3D column vector.');

% Get a euler rotation by pi around the x-axis
[R_pi] = R_x(pi);     
% Rotate the vector through pi
v = R_pi*v;
end