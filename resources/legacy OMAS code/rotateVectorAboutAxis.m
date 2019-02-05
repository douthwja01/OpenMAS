% ROTATE VECTOR THROUGH AN ANGLE, AROUND A GIVEN AXIS
function [newVector] = rotateVectorAboutAxis(oldVector,axisVector,theta)
% This function is designed to calculate a vector
% following a rotation around a given axis vector, through a
% given angle.
% INPUTS:
% oldVector  - The initial vector
% axisVector - The axis of rotation
% theta      - The angle of rotation
% OUTPUTS:
% newVector  - The rotated 3D vector

% INPUT HANDLING
if ~exist('oldVector','var') || ~exist('axisVector','var')
    warning('Both an input vector and rotation axis must be provided.')
    return
end
if length(oldVector) ~= length(axisVector)
    warning('Input vector and axis vector must be the same length.');
    return
end

% NORMALISE THE AXIS VECTOR
axisVector = axisVector/(sqrt(sum(axisVector.^2)));  % Normalize rotation axis

% GET THE CROSS PRODUCT PROECTION BETWEEN THE AXIS AND VECTOR
crossVector = cross(axisVector,oldVector);

% DETERMINE THE MAPPING OF EACH VECTOR COMPONENTS
newVector = cos(theta)*oldVector ...
    + (crossVector)*sin(theta)  ...
    + axisVector*(dot(axisVector,oldVector))*(1 - cos(theta));

end