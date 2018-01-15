%% OPENMAS VISUAL UPDATE CYCLE (OMAS_visuals.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% This utility preforms the visualisation update for OpenMAS as the time
% progression continues.

% Author: James A. Douthwaite 01/06/2017

function [] = OMAS_visuals(METAobject)
% This function attempts to update the current VR scene with the an objects
% META object.

% OBJECT PROPERTIES
position = METAobject.globalState(1:3);
velocity = METAobject.globalState(4:6);
quaternion = METAobject.globalState(7:10);
eulerRotations = simulation_axisTools.quaternionToEuler(quaternion);

% CONVERT TO VR AXIS CONVENTION
positionVR = position;
positionVR(2) = -position(3);
positionVR(3) = position(2);        % The objects global VR position
velocityVR = velocity;
velocityVR(2) = -velocity(3);
velocityVR(3) = velocity(2);        % The objects global VR velocity
rotationVR = zeros(3,1);
rotationVR(1) = -eulerRotations(2);
rotationVR(2) = eulerRotations(3);
rotationVR(3) = -eulerRotations(1); % The VR rotations from the body rotations
% GET THE DCM DESCRIBING THE ATTITUDE UPDATE
DCM_update = simulation_axisTools.eulerToRotationMatrix(rotationVR);

% FROM THE DCM MATRIX, GET THE ROTATED VR-MODEL ROTATION

display('UPDATING VR OBJECT');

% THE OBJECTS VRML MUST BE RELATED TO THE META.OBJECT.objectID

end