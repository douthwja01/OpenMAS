%% GET THE EULER ROTATION MATRIX R(phi,theta,psi) %%%%%%%%%%%%%%%%%%%%%%%%%
function [Reta] = R_eta(phi,theta,psi)

% Input sanity check
if nargin == 3
    eta = [phi,theta,psi]';
elseif nargin == 1
    eta = phi; % phi = [phi,theta,psi]
elseif nargin == 0 
    eta = [sym('phi_t','real');sym('theta_t','real');sym('psi_t','real')];
else
    error('Please provide either a vector of Euler angle or a series of Euler angles.');
end

% Check the inputs
assert(numel(eta) == 3,'Euler angles must be provided by a vector [3x1].');
assert(isa(eta,'sym') || isnumeric(eta),'Euler angles must be a symbolic or numeric value.');

% Get the progressive rotation matrices
R_phi   = R_x(eta(1));
R_theta = R_y(eta(2));
R_psi   = R_z(eta(3));
% Compute the rotations in Z-Y-X order
Reta = R_psi*R_theta*R_phi;
end