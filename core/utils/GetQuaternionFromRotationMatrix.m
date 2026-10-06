% CONVERT ROTATION MATRIX TO QUATERNION
function [q] = GetQuaternionFromRotationMatrix(R)
% This function is designed to convert from a rotation matrix
% to an equivalent quaternion. This function is also parallel
% to "rotm2quat.m".

q = zeros(4,1);
% Assemble the quaternion elements
q(1) = 0.5*sqrt(1 + R(1,1) + R(2,2) + R(3,3));
q(2) = (R(3,2) - R(2,3))/(4*q(1));
q(3) = (R(1,3) - R(3,1))/(4*q(1));
q(4) = (R(2,1) - R(1,2))/(4*q(1));

assert(any(isnan(q)) == 0,'quaternion calculation failed');

% Normalise the quaternion
q = unit(q);
end