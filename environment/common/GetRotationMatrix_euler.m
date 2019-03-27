%% GET THE EULER ROTATION MATRIX R(phi,theta,psi) %%%%%%%%%%%%%%%%%%%%%%%%%
function [R_eta] = GetRotationMatrix_euler(eta)

% Check the inputs
assert(numel(eta) == 3,'Euler angles must be provided by a vector [3x1].');
assert(isa(eta,'sym') || isnumeric(eta),'Euler angles must be a symbolic or numeric value.');
% Get the progressive rotation matrices
R_phi   = GetRotationMatrix_x(eta(1));
R_theta = GetRotationMatrix_y(eta(2));
R_psi   = GetRotationMatrix_z(eta(3));
% Compute the rotations in Z-Y-X order
R_eta = R_psi*R_theta*R_phi;
end