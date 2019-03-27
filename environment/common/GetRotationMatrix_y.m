
function [Ry] = GetRotationMatrix_y(theta)

% Rotate about the y-axis (pitch)
Ry = [cos(theta) 0 sin(theta);
               0 1          0;
     -sin(theta) 0 cos(theta)];
end