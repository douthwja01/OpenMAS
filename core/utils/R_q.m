% GET ROTATION MATRIX FROM QUATERNION
function [R] = R_q(q)
% This function defines the equivalent rotation matrix from the
% provided quaternion. (Associated block "Quaternion to DCM")
% INPUT:
% q    - The quaternion rotation
% OUTPUT:
% R_AB - The rotation matrix through A to B
% R_BA - The reverse rotation matrix from B to A

assert(numel(q) == 4,'Input quaternion has wrong dimensions');

R = zeros(3,3);
% Normalise the quaternion
q = unit(q);                                   

% SAME AS THE ROBOTICS TOOL BOX....
R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;                
R(1,2) = 2*(q(2)*q(3) - q(1)*q(4));
R(1,3) = 2*(q(1)*q(3) + q(2)*q(4));
R(2,1) = 2*(q(2)*q(3) + q(1)*q(4));
R(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
R(2,3) = 2*(q(3)*q(4) - q(1)*q(2));
R(3,1) = 2*(q(2)*q(4) - q(1)*q(3));
R(3,2) = 2*(q(1)*q(2) + q(3)*q(4));
R(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
end