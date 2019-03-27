% CONVERT ROTATION ANGLES INTO A QUATERNION
function [q] = GetQuaternionFromEulers(eulers)
% THis function converts euler angle rotations into a
% quaternion via its components. Associated block:
% "Rotation Angles to Quaternions".
eulers = 0.5*eulers;
% Build vector of trigonometric arguements
trigArgs = [sin(eulers);cos(eulers)];
% Assemble Quaternion components
q = zeros(4,1);
q(1) = trigArgs(4)*trigArgs(5)*trigArgs(6) + trigArgs(1)*trigArgs(2)*trigArgs(3);
q(2) = trigArgs(4)*trigArgs(5)*trigArgs(3) - trigArgs(1)*trigArgs(2)*trigArgs(6);
q(3) = trigArgs(4)*trigArgs(2)*trigArgs(6) + trigArgs(1)*trigArgs(5)*trigArgs(3);
q(4) = trigArgs(1)*trigArgs(5)*trigArgs(6) - trigArgs(4)*trigArgs(2)*trigArgs(3);
% Renormalise
q = unit(q);
end