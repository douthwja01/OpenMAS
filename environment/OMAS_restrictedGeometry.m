%% RESTICTED geometry DEFINITION (OMAS_restrictedGeometry.m) %%%%%%%%%%%%%%
% This function computes the geometry that can be observed within a defined
% radius.

function [restrictedGeometry] = OMAS_restrictedGeometry(point,radius,geometry)
% INPUTS:
% centroid - The 3D cartesian reference point.
% radius   - The detection radius of the observational sphere.
% geometry - The geometry of the second object in the reference frame of
%            the centroid.

% TO ENSURE FUNCTIONALITY
sendCompleteGeometryOnDetection = 1;

newSurfaceMembers = [];
% EVALUATE THE geometry'S FACE MEMBERS
for face = 1:size(geometry.faces,1)
    % NOTE:
    % - If all vertices are outside of the sphere, the sphere may still 
    %   intersect the plane at an edge or mid-plane.
    
    % MEMBERS IF THE PLANE
    memberID_A = geometry.faces(face,1);
    memberID_B = geometry.faces(face,2);
    memberID_C = geometry.faces(face,3);
    % TEST THE FACE AGAINST THE SPHERE
    faceIsInsideRadius = OMAS_sphereTriangleIntersection(point,radius,...
                                                         geometry.vertices(memberID_A,:)',...
                                                         geometry.vertices(memberID_B,:)',...
                                                         geometry.vertices(memberID_C,:)');
    if sendCompleteGeometryOnDetection && faceIsInsideRadius
    % IF ANY FACE IS VISIBLE, SEND THE COMPLETE GEOMETRY    
        restrictedGeometry = geometry;
        return
    elseif faceIsInsideRadius
    % IF THE FACE IS FOUND TO VIOLATE THE CONSTRAINT
%         fprintf('\t face %d interacts with the sphere\n',face);
        % add all vertices belonging to that face to the matrix
        newSurfaceMembers = vertcat(newSurfaceMembers,memberID_A,memberID_B,memberID_C);
    end
end
% PARSE COMMON POINTS
newSurfaceMembers = unique(newSurfaceMembers,'rows');   % The unique vertices
newVertices = geometry.vertices(newSurfaceMembers,:);   % Extract the common vertices

% THE CONSTRAINT SURFACE
% sphereGeometry = OMAS_graphics.defineSphere(zeros(3,1),radius,6);
% [A B] = SurfaceIntersection(sphereGeometry,geometry);


% TRIANGULATE THE NEW VOLUME
triangulatedVolume = triangulation(boundary(newVertices(:,1),newVertices(:,2),newVertices(:,3),1),newVertices);
restrictedGeometry = geometry;
restrictedGeometry.vertices = triangulatedVolume.Points;
restrictedGeometry.faces = triangulatedVolume.ConnectivityList;
restrictedGeometry.normals = OMAS_graphics.normals(restrictedGeometry);    % Compute the surface normals

% DEBUG
% fig = figure(2);
% ax = axes(fig);
% axis square vis3d;
% % THE OBJECT GEOMETRY
% patch(ax,...
%     'Faces',geometry.faces,...
%     'Vertices',geometry.vertices,...
%     'faceColor','g');
% % THE CONSTRAINT VOLUME
% patch(ax,...
%     'Faces',sphereGeometry.faces,...
%     'Vertices',sphereGeometry.vertices,...
%     'faceColor','b',...
%     'faceAlpha',0.2);
% % THE RESULTING GEOMETRY
% patch(ax,...
%     'Faces',restrictedGeometry.faces,...
%     'Vertices',restrictedGeometry.vertices,...
%     'faceColor','b',...
%     'faceAlpha',0.2);
% end

end

% THE SQUARED DISTANCE TO ALL VERTICES
% verticesSq = (geometry.vertices(:,1) - point(1)).^2 ...
%            + (geometry.vertices(:,2) - point(2)).^2 ...
%            + (geometry.vertices(:,3) - point(3)).^2;
% % VERTEX INSIDE RADIUS LOGICALS
% vertexDetectionLogicals = verticesSq <= radius^2;
% vertexIDset = 1:1:numel(vertexDetectionLogicals);
% 
% % VISIBLE POINT ID's
% visibleIDs = vertexIDset(vertexDetectionLogicals);
% if numel(visibleIDs) < 1
%     disp('no vertices are visible');
% end

% [A,~] = ismember(relativeFaces,vertIndices);
% completeFaceIndices = sum(A,2) > 2;                        % Minimum number of connections

