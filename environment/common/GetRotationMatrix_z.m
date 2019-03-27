
function [Rz] = GetRotationMatrix_z(psi)

% Rotate about the z-axis (yaw)
Rz = [cos(psi) -sin(psi) 0;
      sin(psi)  cos(psi) 0;
             0         0 1];
end