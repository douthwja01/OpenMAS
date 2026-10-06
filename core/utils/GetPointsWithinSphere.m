function [vertexIndices,distSq] = GetPointsWithinSphere(point,radius,xyzData)

% CHECK INPUTS
assert(numel(point) == 3,   'Point must be specified as a 3D point [3x1]');
assert(size(xyzData,2) == 3,'Vertex data must be specified a list of 3D vertices [nx3]');
assert(numel(radius) == 1,  'The radius is a scalar comparitor.');

% DEFAULT VALUES
% internalPoints = xyzData;               % Return the full set
vertexIndices = 1:1:size(xyzData,1);      % Return a substitute index vector
distSq = inf(1,3);

% DEFAULT RESPONSE
if isinf(radius)
    return
end
% VECTORISED SQUARED DISTANCE CALCULATION
distSq = (xyzData(:,1) - point(1)).^2 + (xyzData(:,2) - point(2)).^2 + (xyzData(:,3) - point(3)).^2;
% SQUARED DISTANCE COMPARITOR
ind = distSq < radius^2;
% THE REDUCED POINT CLOUD INSIDE THE RADIUS
[vertexIndices,~] = find(ind); 
end