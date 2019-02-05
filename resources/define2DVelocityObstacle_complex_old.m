% DECOMPOSE THE OBSTACLE GEOMETRIES INTO A VO SET
function [VO] = define2DVelocityObstacle_complex(obj,p_a,v_a,r_a,p_b,v_b,geometry,tau)
% This function computes the VO set from a list of complex
% obstacles defined by their .GEOMETRY parameter.

% NOTES:
% - The .geometry parameter is the geometry of the object pre-
%   positioned relative to the agent, rotated and scaled.

% INPUT HANDLING
if ~obj.VIRTUAL.is3D
    p_a = [p_a;0]; v_a = [v_a;0];
    p_b = [p_b;0]; v_b = [v_b;0];
end
% CONSTANTS
origin   = zeros(3,1);
VOorigin = v_b;
lambda   = p_b - p_a;

% ////////////////// PROBLEM PARAMETERS ///////////////////////
planeNormal = [0;0;1];
unit_reference = [1;0;0];
lambda_planar = OMAS_geometry.vectorPlanarProjection(planeNormal,lambda); % planar projection of position
norm_lambda_planar = norm(lambda_planar);
unit_lambda_planar = lambda_planar/norm_lambda_planar;
% PREPARE OBJECT VERTICES
lambda_vertices = geometry.vertices - origin';  % The relative vertex coordinates
lambda_vertices(:,3) = 0;
lambda_vertices = unique(lambda_vertices,'rows');
% THE DIRECTION OF ROTATION FROM THE FORWARD VECTOR
%             centroidIsRightLogical = -sign(dot(cross(unit_reference,unit_lambda_planar),planeNormal)); % 1 (right) -1 (left) % Direction of the rotation relative to the centerline
hold on;
% CENTER-LINE
q = quiver3(gca,origin(1),origin(2),origin(3),lambda_planar(1),lambda_planar(2),lambda_planar(3),'-.');
q.AutoScaleFactor = 1;
% /////////////////////////////////////////////////////////////
for i = 1:size(lambda_vertices,1)
    % GET THE VERTEX PARAMETERS
    vertex_planar = lambda_vertices(i,:)';
    norm_vertex = norm(vertex_planar);
    unit_vertex = vertex_planar/norm_vertex;
    
    % VERTEX IS ON RIGHT/LEFT OF THE FORWARD VECTOR
    vertexIsOnRightLogical = -sign(dot(cross(unit_reference,unit_vertex),planeNormal));       % Direction of the rotation relative to the centerline
    
    % THE MINIMUM VECTOR FROM THE CURRENT VERTEX
    %                 perpVector = r_a*cross(planeNormal,unit_vertex);
    
    % /////////////////////////////////////////////////////////
    
    % THE VERTEX MODIFIER (r_a) ANGLE
    %                 alphaAngles(i) = asin(r_a/norm_vertex);
    vertexAngle = -vertexIsOnRightLogical*acos(dot(unit_reference,unit_vertex));
    if abs(vertexAngle) == pi
        %                     vertexAngle = 0;
        continue
    end
    % BUILD VERTEX ANALYSIS STRUCTURE
    vertexData(i) = struct('point',vertex_planar,...
        'distance',norm_vertex,...
        'perp',r_a*cross(planeNormal,unit_vertex),...
        'isOnRight',vertexIsOnRightLogical,...
        'angle',vertexAngle);
    
    % //////////////////// DEBUG PLOTS ////////////////////////
    if vertexData(i).isOnRight > 0
        scatter3(gca,vertex_planar(1),vertex_planar(2),vertex_planar(3),'ro');
    else
        scatter3(gca,vertex_planar(1),vertex_planar(2),vertex_planar(3),'go');
    end
end
% All the vertex data is now held in the structure 'vertexData'
candidateAngles = [vertexData(:).angle]
[~,rightInd]  = min(candidateAngles);
rightVertex = vertexData(rightInd);
[~,leftInd] = max(candidateAngles);
leftVertex = vertexData(leftInd);


% CALCULATE THE VO PARAMETERS
leftEdgeAngle  = leftVertex.angle   %- %asin(r_a/leftVertex.distance);
rightEdgeAngle = rightVertex.angle  %+ %asin(r_a/rightVertex.distance);

% AUGMENT BY THE MINIMUM SEPARATION RADIUS
%             leftVertex.point + leftVertex.perp


equivalentOpenAngle  = rightEdgeAngle - leftEdgeAngle;
unit_leadingTangent  = OMAS_geometry.rodriguesRotation(unit_reference,planeNormal,rightEdgeAngle);
unit_trailingTangent = OMAS_geometry.rodriguesRotation(unit_reference,planeNormal,leftEdgeAngle);

% /////////////// VO IS DEFINED, COMPUTE LOGIC ////////////////

% DETERMINE WHERE Va LIES
if dot(unit_reference,unit_leadingTangent) > dot(unit_reference,unit_trailingTangent)
    isVaLeading = 1;
else
    isVaLeading = 0;
end
% DETERMINE IF Va IS INSIDE THIS VO
isVaInsideCone = 0;
VOtolerance = 1E-8;

candidateVector = unit_reference - VOorigin;

% EVALUATE THE PROJECTION AGAINST LEADING EDGE FIRST
VOprojection   = norm(candidateVector)*cos(abs(leftEdgeAngle));
candProjection = unit_lambda_planar'*candidateVector;
projDiff = (candProjection - VOprojection);
if projDiff > VOtolerance && isVaLeading
    isVaInsideCone = 1;
end
% EVALUATE THE PROJECTION AGAINST TRAILING EDGE
VOprojection   = norm(candidateVector)*cos(abs(rightEdgeAngle));
candProjection = unit_lambda_planar'*candidateVector;
projDiff = (candProjection - VOprojection);
if projDiff > VOtolerance && ~isVaLeading
    isVaInsideCone = 1;
end

effectiveRadii = 5;

% DEFINE THE VO STRUCTURE
VO.apex                     = origin; % VOorigin;
VO.axisUnit                 = unit_lambda_planar;
VO.axisLength               = norm_lambda_planar;
VO.openAngle                = equivalentOpenAngle;
VO.leadingEdgeUnit          = unit_leadingTangent;
VO.trailingEdgeUnit         = unit_trailingTangent;
VO.isVaLeading              = isVaLeading;
VO.isVaInsideCone           = isVaInsideCone;
VO.truncationTau            = tau;
VO.truncationCircleCenter   = (lambda_planar)/tau + VOorigin;
VO.truncationCircleRadius   = (effectiveRadii + r_a)/tau;

% DEBUG PLOT
CCorigin = zeros(3,1);
scale = 10;
q.AutoScaleFactor = 1;
q = quiver3(gca,CCorigin(1),CCorigin(2),CCorigin(3),...
    VO.leadingEdgeUnit(1)*scale,VO.leadingEdgeUnit(2)*scale,VO.leadingEdgeUnit(3)*scale,'g');
q.AutoScaleFactor = 1;
q = quiver3(gca,CCorigin(1),CCorigin(2),CCorigin(3),...
    VO.trailingEdgeUnit(1)*scale,VO.trailingEdgeUnit(2)*scale,VO.trailingEdgeUnit(3)*scale,'r');
q.AutoScaleFactor = 1;

% CONVERT TO 2D
VO.apex = VO.apex(1:2,1);
VO.axisUnit = OMAS_geometry.unit(VO.axisUnit(1:2,1));
VO.leadingEdgeUnit  = OMAS_geometry.unit(VO.leadingEdgeUnit(1:2));
VO.trailingEdgeUnit = OMAS_geometry.unit(VO.trailingEdgeUnit(1:2));
VO.truncationCircleCenter = VO.truncationCircleCenter(1:2);
end