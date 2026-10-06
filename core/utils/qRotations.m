% CONVERT QUATERNION INTO ROTATION ANGLES
function [eulerAngles] = qRotations(q)                    
% Gets the rotation angles from an equivalent quaternion
% attitude. Associated block:
% "Quaternions to Rotation Angles"
assert(size(q,1) == 4,'A quaternion column vector is expected 4x1');

% Normalise the quaternion
q = unit(q);
eulerAngles = zeros(3,1);
% QUATERNION NOTATION B [roll pitch yaw]
eulerAngles(1) = atan2(2*(q(1)*q(2) + q(3)*q(4)),(1 - 2*(q(2)^2 + q(3)^2)));              % Phi   (roll)
eulerAngles(2) = asin( 2*(q(1)*q(3) - q(4)*q(2)));             % Theta (pitch)
eulerAngles(3) = atan2(2*(q(1)*q(4) + q(2)*q(3)),(1 - 2*(q(3)^2 + q(4)^2)));              % Psi   (yaw)
end