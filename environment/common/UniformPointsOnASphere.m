
function [points] = UniformPointsOnASphere(p,r)
% This function generates a point cloud equally distributed across the
% surface of a sphere.

% Input sanity check
if nargin < 2
    r = 1;
end

% Constants
viewDir = 3;
goldenRatio     = (1 + sqrt(5))/2;
angleIncrement  = pi*2*goldenRatio;

points = [];
for i = 1:numel(viewDir)
    t = i/viewDir;
    inclination = acos(1-2*t);
    azimuth = angleIncrement*i;
    
    % Convert spherical to cartesian
    x = p(1) + r*sin(inclination)*cos(azimuth);
    y = p(2) + r*sin(inclination)*sin(azimuth);
    z = p(3) + r*cos(inclination);
    
    % Assign point to matrix
    points(i,:) = [x,y,z];
end


end

%     static BoidHelper () {
%         directions = new Vector3[BoidHelper.numViewDirections];
% 
%         float goldenRatio = (1 + Mathf.Sqrt (5)) / 2;
%         float angleIncrement = Mathf.PI * 2 * goldenRatio;
% 
%         for (int i = 0; i < numViewDirections; i++) {
%             float t = (float) i / numViewDirections;
%             float inclination = Mathf.Acos (1 - 2 * t);
%             float azimuth = angleIncrement * i;
% 
%             float x = Mathf.Sin (inclination) * Mathf.Cos (azimuth);
%             float y = Mathf.Sin (inclination) * Mathf.Sin (azimuth);
%             float z = Mathf.Cos (inclination);
%             directions[i] = new Vector3 (x, y, z);
%         }
%     }
    
    