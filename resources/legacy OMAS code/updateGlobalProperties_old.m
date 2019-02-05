% UPDATE THE GLOBAL REPRESENTATION (WITH WAYPOINT COMPLETION)
function [obj] = updateGlobalProperties(obj,dt,localState_update)
% This function is used to calculate the new global (.VIRTUAL)
% parameters for the current object.
% INPUTS:
% obj   - The current object class
% dt    - The time since this function was last ran
% localState_update - The objects new state
% OUTPUTS:
% .globalPosition - The new global position
% .globalVelocity - The new global velocity
% .quaternion     - The new quaternion pose
% .R              - The new rotation matrix

% updates the agents global properties from
% X_global(@t=k] to X_global(@t=k+1)

% CONSTANT PARAMETERS
rotationIndices = 4:6;
velocityIndices = 7:9;

% IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
    localState_update(7:end) = NaN([6,1]);
    obj.VIRTUAL.idleStatus = 1;
end

% DEFINE UPDATE PARAMETERS
localState_FRD_k_plus = localState_update;                     % Collect parameters for this time step
globalPosition_k = obj.VIRTUAL.globalPosition;
globalVelocity_k = obj.VIRTUAL.globalVelocity;

% //////// GET THE UPDATED QUATERNION /////////////////////////
[localState_FLU_k_plus] = OMAS_axisTools.convertNEDstateToENU(localState_FRD_k_plus);
[localState_XYZ_k_plus] = OMAS_axisTools.convertENUStateToMatlab(localState_FLU_k_plus);

% EXTRACT ROTATIONS
XYZrotations_k_plus = localState_XYZ_k_plus(rotationIndices);

% REVERSE ROTATION ORDER
ZYXrotations_k_plus = transpose([0,0,1;0,1,0;1,0,0]*XYZrotations_k_plus);
[quaternion_k_plus] = eul2quat(ZYXrotations_k_plus,'ZYX');
quaternion_k_plus = transpose(quaternion_k_plus);
% NEW ROTATION MATRIX
[R_k_plus] = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
% ////////// UPDATE GLOBAL POSE ///////////////////////////////
globalVelocity_k_plus = R_k_plus*localState_XYZ_k_plus(velocityIndices);
globalPosition_k_plus = globalPosition_k + dt*globalVelocity_k;
% ////////// REASSIGN K+1 PARAMETERS //////////////////////////
obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
obj.VIRTUAL.R = R_k_plus;
obj.localState = localState_FRD_k_plus;                        % Reassign the obj.localstate
end