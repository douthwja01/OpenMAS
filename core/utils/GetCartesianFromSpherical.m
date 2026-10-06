% CONVERT SPHERICAL TO CARTESIAN
function [p] = GetCartesianFromSpherical(r,psi,theta)
% This function calculates the relative position from a spherical
% coordinate system that assumes angles are measured from the
% agents heading vector Xref.
% INPUTS:
% range     - The objects radial seperation (m)
% azimuth   - The objects relative azimuth angle (rad)
% elevation - The objects relative elevation (rad)
% OUTPUTS:
% cartesianPosition  - The new 3D position in local coordinates (m)

% DEFINE THE POSITION AS A VECTOR INTERVAL
p = [cos(psi)*cos(theta);...
     sin(psi)*cos(theta);...
              sin(theta)]*r;
end