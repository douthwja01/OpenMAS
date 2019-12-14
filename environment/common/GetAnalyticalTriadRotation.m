% ANALYTICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
function [q] = GetAnalyticalTriadRotation(referenceTriad,targetTriad)
% This function defines the quaternion describing the rotations
% between a reference triad and a second triad.

% NORMALISE THE TRIAD VECTORS
for dim = 1:size(referenceTriad,2)
    [referenceTriad(:,dim)] = unit(referenceTriad(:,dim));
    [targetTriad(:,dim)] = unit(targetTriad(:,dim));
end

% EXTRACT REFERENCE AXES
xAxis = targetTriad(:,1);
zAxis = targetTriad(:,3);   % Get the rotated body axes
xAxis_ref = referenceTriad(:,1);
zAxis_ref = referenceTriad(:,3); % Get the reference unit ENU triad (trying to get to)

% FIRST ALIGN THE Z-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE QUATERNION ROTATION TO ALLIGN THE Z-AXES
[q_zAlign] = qArgument(zAxis,zAxis_ref);                                   % Quaternion aligning global and body z vectors
% [R_zAlign] = OMAS_geometry.quaternionToRotationMatrix(q_zAlign);         % Equivalent rotation matrix
[R_zAlign] = R_q(q_zAlign);                                % Equivalent rotation matrix
% ALIGN THE X AXIS IN THE Z-PLANE
xAxis_intermediate = R_zAlign*xAxis;
% TAKE ITS PROJECTIONS IN THE XY PLANE & RENORMALISE
xAxis_intermediate(3) = 0;
xAxis_intermediate = unit(xAxis_intermediate);
% OTHERWISE JUST ALIGN THE X-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE QUATERNION ROTATION TO ALLIGN THE X-AXES
[q_xAlign] = qArgument(xAxis_intermediate,xAxis_ref);
% [R_xAlign] = OMAS_geometry.quaternionToRotationMatrix(q_xAlign);
[R_xAlign] = R_q(q_xAlign);

% COMPUTE THE COMPOSITE ROTATION MATRIX
comp_rotation = R_xAlign * R_zAlign;
% COVERT THE ROTATION MATRIX TO A QUATERNION DESCRIBING THE
% ROTATION FROM TRIAD_REF TO TRIAD_FINAL
q = rotm2quat(comp_rotation)';

% RENORMALISE
q = unit(q);
% OUTPUT CATCHA
assert(any(isnan(q)) == 0,'No valid quaternion found');
end