% ROTATE A VECTOR THROUGH A QUATERNION
function [v1] = qVectorRotation(q,v0)
% This function rotates a cartesian vector through a quaternion
% rotation. Associated block:
% "Quaternion Rotation"
if numel(q) ~= 4
    error(char(['Input quaternion is of length %d ',length(q)']));
end
% Normalise the quaternion rotation
q = unit(q);
% Rotate the vector through the quaternion elements
v1 = 2*[(0.5 - q(3)^2 - q(4)^2), (q(1)*q(4) + q(2)*q(3)), (q(2)*q(4) - q(1)*q(3));...
        (q(2)*q(3) - q(1)*q(4)), (0.5 - q(2)^2 - q(4)^2), (q(1)*q(2) + q(3)*q(4));...
        (q(1)*q(3) + q(2)*q(4)), (q(3)*q(4) - q(1)*q(2)), (0.5 - q(2)^2 - q(3)^2)]*v0;
end