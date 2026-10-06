% GET THE QUATERNION BETWEEN TWO VECTORS
function [q] = qArgument(u,v)
q = zeros(4,1);
% Normalise the quaternion
u = u/norm(u);
v = v/norm(v);
% Get the axis vector
q(2:4) = cross(u,v);
% Define the rotation about that vector
q(1) = sqrt((norm(u)^2)*(norm(v)^2)) + dot(u,v);
% Normalise the quaternion
[q] = unit(q);
end