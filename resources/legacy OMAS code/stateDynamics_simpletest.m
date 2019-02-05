% UPDATE AGENT STATE VECTOR FROM A VELOCITY VECTOR
function dXdt = stateDynamics_simple(X,velocity_k_plus,omega_k_plus)
% This function assumes that the velocity changes are
% implemented this timestep directly, integration then occurs
% including the updates from this timestep:
% State vector is assumed:
% 2D - [x y psi dx dy dpsi]'
% 3D - [x y z phi theta psi dx dy dz dphi dtheta dpsi]'

% CHECK EXPECTED INPUT DIMENSIONS
Aflag = numel(velocity_k_plus) == 3 && numel(omega_k_plus) == 3; % IS 3D
Bflag = numel(velocity_k_plus) == 2 && numel(omega_k_plus) == 1; % IS 2D
assert((Aflag + Bflag) == 1,'Integrator inputs are of the wrong dimensions');

% Either rotation occurs in one dimension (2D) or in three (3D)
rotationDOF = numel(omega_k_plus); % Is heading a yaw (2D) or [phi,theta psi] (3D)
if rotationDOF == 1
    linearIndices = 1:2;
    angularIndices = 3;
else
    linearIndices = 1:3;
    angularIndices = 4:6;
end

% CALCULATE THE STATE UPDATE
% X = zeros(2*stateNumber,1);
% X(linearIndices)  = X(linearIndices)  + dt*velocity_k_plus;
% X(angularIndices) = X(angularIndices) + dt*omega_k_plus;
% X(stateNumber+linearIndices)  = velocity_k_plus;
% X(stateNumber+angularIndices) = omega_k_plus; % Record the input velocities

stateNum = numel([linearIndices,angularIndices]);

% THE STATE DIFFERENTIAL
dXdt = zeros(2*stateNum,1);
dXdt(linearIndices) = velocity_k_plus;
dXdt(angularIndices) = omega_k_plus;
dXdt(stateNum+linearIndices) = velocity_k_plus - X(stateNum+linearIndices); % Infer the accelerations
dXdt(StateNum+angularIndices) = omega_k_plus - X(StateNum+angularIndices);
end