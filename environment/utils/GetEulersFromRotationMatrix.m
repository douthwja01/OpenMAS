% EULER ROTATION MATRIX TO EULER ANGLES
function [eulers] = GetEulersFromRotationMatrix(R)
% Rotation of the fixed to the global pose
eulers = zeros(3,1);
eulers(1) =  atan2(R(3,2),R(3,3));  % The roll angle
eulers(2) = -asin(R(3,1));             % The pitch angle
eulers(3) =  atan2(R(2,1),R(1,1));  % The yaw angle
end