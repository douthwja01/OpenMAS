
function [Rx] = GetRotationMatrix_x(phi)

% Rotate about the x-axis (roll)
Rx  = [ 1        0         0;
        0 cos(phi) -sin(phi);
        0 sin(phi)  cos(phi)];
end