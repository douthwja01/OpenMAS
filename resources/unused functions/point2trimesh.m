%% POINT TO TRI-MESH (point2trimesh.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the point 



function [ distances, surface_points, faces2, vertices2, corresponding_vertices_ID, new_faces_ID ] = point2trimesh( varargin )
% -------------------------------------------------------
%
%    point2trimesh - Distance between a point and a triangulated surface in 3D 
%    
%    The shortest line connecting a point and a triangulation in 3D is 
%    computed. The nearest point on the surface as well as the distance
%    is returned. The distance is signed according to face normals to
%    identify on which side of the surface the query point resides. 
%    The implementation is optimized for speed, and depending on your
%    application you can use linear or parallel computation. 
%
%    Point insertion functionality 
%    (this feature is experimental and not optimized for speed):
%    If the function is called with more than two output arguments, 
%    the surface points are included into the given triangulation 
%    and Delaunay conditions are restored locally. If triangles with small
%    angles occur, additional vertices are inserted to eliminate them
%    if possible. 
%
%    Algorithm: From every query point, 
%       - the nearest vertex
%       - the nearest point on the edges and
%       - the nearest point on the triangle's surfaces
%    is calculated and the minimum distance out of these three is returned.   
%
%    Ver. 1.0.0
%
%    Created:         Daniel Frisch        (29.03.2015)
%    Last modified:   Daniel Frisch        (24.06.2015)
%
% ------------------------------------------------------
%
%  Inputs:
%      Pairs of parameter names and corresponding values. 
%      Structure arrays are expanded into separate inputs, 
%      where each field name corresponds to an input parameter name. 
%      Parameters:
%      - 'Faces'         (#faces    x 3) Triangulation connectivity list. Each row defines a face, elements are vertex IDs.  
%      - 'Vertices'      (#vertices x 3) Point matrix. Columns are x,y,z coordinates, row numbers are vertex IDs.  
%      - 'QueryPoints'   (#qPoints  x 3) Columns are x,y,z coordinates; each row defines a query point. Can be empty. 
%      - 'MaxDistance'   (1 x 1)         If the distance between a surface_point and its nearest vertex is within this range, 
%                                        no new vertex is inserted into the mesh. This helps avoiding
%                                        triangles with small angles. (default: 1/10 the smallest inradius)
%      - 'UseSubSurface' (1 x 1)         Logical. If true (default), the distance to edges and surfaces is only calculated
%                                        for faces that are connected to the vertex nearest to the query point.  
%                                        This speeds up the calculation but if the distance between two opposite parts 
%                                        of the surface is less than the spacing of the vertices, wrong results are produced.
%                                        In the vectorized algorithm, 'SubSurface' is always false.  
%      - 'Algorithm'     'linear' (default): query points are processed successively in a 'for' loop. 
%                            Use this if you have only few query points. 
%                        'parallel': query points are processed in parallel with 'parfor'. 
%                            If no parallel pool exists, Matlab creates one automatically (which takes half a minute).
%                            It shuts down after 30 min idle time by default, but you can change that in the "Parallel Preferences". 
%                            If Matlab doesn't correctly recognize the number of cores of your processor, 
%                            change the "Number of workers" in "Manage Cluster Profiles".   
%                         'vectorized': query points are processed altogether in a vectorized manner.
%                            Be careful, this needs much RAM if you have many query points. 
%                            This option is mostly included to show that vectorization does not always speed up the calculation:
%                            In the linear code, every query point can be assigned an individual cutout of the whole surface. 
%                         'linear_vectorized_subfunctions': query points and their individual cutout surfaces
%                            are processed successively in a for loop by functions that are capable of processing more than one point. 
%                            This option is included to show that the non-vectorized functions are faster 
%                            if only one point is processed at a time. 
%                         'parallel_vectorized_subfunctions': query points and their individual cutout surfaces
%                            are processed in parallel in a parfor loop by functions that are capable of processing more than one point. 
%                            Again, this option is included to show that the non-vectorized functions are faster 
%                            if only one point is processed at a time. 
%
%  Outputs:
%      - distances      (#qPoints   x 1)   Vector with the point-surface distances; sign depends on normal vectors. 
%      - surface_points (#qPoints   x 3)   Matrix with the corresponding nearest points on the surface. 
%      - faces2         (#faces2    x 3)   Connectivity matrix of the triangulation including the surface_points as dedicated vertices  
%      - vertices2      (#vertices2 x 3)   Point/Vertex matrix of the triangulation including the surface_points as dedicated vertices 
%      - corresponding_vertices_ID    (#qPoints x 1) Vector with the IDs of the vertices corresponding to the query points
%      - new_faces_ID                 Vector with the IDs of the new or modified faces (to give them a different color, for example)
%
%
% Usage example:
%      FV.faces    = [5 3 1; 3 2 1; 3 4 2; 4 6 2];
%      FV.vertices = [2.5 8.0 1; 6.5 8.0 2; 2.5 5.0 1; 6.5 5.0 0; 1.0 6.5 1; 8.0 6.5 1.5];
%      points      = [2 4 2; 4 6 2; 5 6 2];
%      [distances,surface_points] = point2trimesh(FV, 'QueryPoints', points); 
%      patch(FV,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
%      plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
%      plot3M(points,'*r')
%      plot3M(surface_points,'*k')
%      plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);shiftdim(points,-1)*NaN],[],3),'k')
%
%


%% Parse Inputs 

% valdiation functions
faceChk   = @(x) validateattributes(x,{'double','int32' },{'real','finite','nonnan','positive'   ,'integer','size',[NaN 3]});
vertChk   = @(x) validateattributes(x,{'double','single'},{'real','finite','nonnan'                        ,'size',[NaN 3]});
pointChk  = @(x) validateattributes(x,{'double'         },{'real','finite','nonnan'                                       });
distChk   = @(x) validateattributes(x,{'double'         },{'real','finite','nonnan','nonnegative'          ,'scalar'      });
logicChk  = @(x) validateattributes(x,{'logical'        },{'scalar'});
charChk   = @(x) validateattributes(x,{'char'           },{'nonempty','vector'});

parser = inputParser;
parser.FunctionName = mfilename;
parser.addParameter('Faces'         ,[]      ,faceChk);
parser.addParameter('Vertices'      ,[]      ,vertChk);
parser.addParameter('QueryPoints'   ,[]      ,pointChk);
parser.addParameter('MaxDistance'   ,[]      ,distChk);
parser.addParameter('UseSubSurface' ,true    ,logicChk);
parser.addParameter('Algorithm'     ,'linear',charChk);

parser.parse(varargin{:});

faces    = double(parser.Results.Faces);
vertices = double(parser.Results.Vertices);
assert(~isempty(faces) && ~isempty(vertices), 'Invalid argument: ''Faces'' and ''Vertices'' mustn''t be empty.')
assert(max(faces(:))<=size(vertices,1), 'The value of ''Faces'' is invalid: the maximum vertex ID is bigger than the number of vertices in ''Vertices''')

qPoints = parser.Results.QueryPoints;
useSubSurface = parser.Results.UseSubSurface;
if nargout>2, insertPoints = true;
else          insertPoints = false; end

% Calculate normals
r1 = vertices(faces(:,1),:);  % (#faces x 3) % 1st vertex of every face
r2 = vertices(faces(:,2),:);  % (#faces x 3) % 2nd vertex of every face
r3 = vertices(faces(:,3),:);  % (#faces x 3) % 3rd vertex of every face
normals = cross((r2-r1),(r3-r1),2); % (#faces x 3) normal vector of every face
normals = bsxfun(@rdivide,normals,sqrt(sum(normals.^2,2))); % (#faces x 3) normalized normal vector of every face

if isempty(qPoints)
    distances = [];
    surface_points = [];
    faces2 = faces;
    vertices2 = vertices;
    corresponding_vertices_ID = [];
    new_faces_ID = [];
    return
end


%% Distance Calculation

nQPoints = size(qPoints,1);

D = NaN(nQPoints,1);
P = NaN(nQPoints,3);

if insertPoints
    max_distance = parser.Results.MaxDistance;
    if isempty(max_distance)
        tri = triangulation(faces,vertices);
        [~,r] = tri.incenter;
        max_distance = min(r)/10;
    end
    max_distance = max(max_distance,100*eps);
    is_new_face = zeros(size(faces,1),1);
    is_new_vertex = zeros(size(vertices,1),1);
    corresponding_vertices_ID = NaN(size(qPoints,1),1);
end

switch parser.Results.Algorithm
    case {'linear','normal'}
        for r = 1:nQPoints
            % Determine the surface points
            point = qPoints(r,:); % (1 x 3) query point
            [d,p,f] = processPoint(faces,vertices,point,normals, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);
            D(r) = d;
            P(r,:) = p;
            if insertPoints
                % Include the surface points as new vertices into the mesh and restore Delaunay conditions 
                [ tri, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( triangulation(faces,vertices), p, f, max_distance, is_new_face, is_new_vertex );
                faces    = tri.ConnectivityList;
                vertices = tri.Points;  
                normals  = tri.faceNormal;
                corresponding_vertices_ID(r) = new_vertex_ID;
            end
        end
        
    case 'parallel'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        parfor r = 1:nQPoints
            point = qPoints(r,:); % (1 x 3) query point
            [d,p] = processPoint(faces,vertices,point,normals, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);
            D(r) = d;
            P(r,:) = p;
        end
        
    case 'vectorized'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        if useSubSurface && ~ismember('UseSubSurface',parser.UsingDefaults)
            warning('You specified ''UseSubSurface'',true, but ''Algorithm'',''vectorized'' always searches on the complete surface') 
        end
        [D1,P1] = distance_to_vertices_vectorized(faces,vertices,qPoints,normals); % (#qPoints x 1), (#qPoints x 3), (#qPoints x 1) 
        [D2,P2] = distance_to_edges_vectorized   (faces,vertices,qPoints,normals);
        [D3,P3] = distance_to_surfaces_vectorized(faces,vertices,qPoints,normals);
        % find minimum distance type
        D = [D1,D2,D3];      % (#qPoints x 3)
        P = cat(3,P1,P2,P3); % (#qPoints x xyz x 3)
        [~,I] = min(abs(D),[],2);
        D = D(sub2ind(size(D),(1:length(I))',I));
        % extract nearest point on surface
        P = permute(P,[2,3,1]); % (xyz x 3 x #qPoints)
        sz = [size(P) 1 1];
        P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (xyz x #qPoints)
        P = P'; % (#qPoints x xyz)
        
    case 'linear_vectorized_subfunctions'
        for r = 1:nQPoints
            % Determine the surface points
            point = qPoints(r,:); % (1 x 3) query point
            [d,p,f] = processPoint(faces,vertices,point,normals, @distance_to_vertices_vectorized,@distance_to_edges_vectorized,@distance_to_surfaces_vectorized, useSubSurface);
            D(r) = d;
            P(r,:) = p;
            if insertPoints
                % Include the surface points as new vertices into the mesh and restore Delaunay conditions 
                [ tri, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( triangulation(faces,vertices), p, f, max_distance, is_new_face, is_new_vertex );
                faces    = tri.ConnectivityList;
                vertices = tri.Points;  
                normals  = tri.faceNormal;
                corresponding_vertices_ID(r) = new_vertex_ID;
            end
        end
        
    case 'parallel_vectorized_subfunctions'
        assert(~insertPoints,'''Algorithm'', ''%s'' doesn''t support including the surface points into the geometry. Call point2trimesh with fewer output arguments or use ''linear'' algorithm.',parser.Results.Algorithm)
        parfor r = 1:nQPoints
            point = qPoints(r,:); % (1 x 3) query point
            [d,p] = processPoint(faces,vertices,point,normals, @distance_to_vertices_vectorized,@distance_to_edges_vectorized,@distance_to_surfaces_vectorized, useSubSurface);
            D(r) = d;
            P(r,:) = p;
        end
        
    otherwise
        error('The value of ''Algorithm'' is invalid.')
end


if insertPoints
    % Despite the Delaunay condition is fulfilled, the
    % triangulation might still contain triangles with small angles.
    % To prevent this, insert vertices at the circumcenters  of these triangles.
    tri = triangulation(faces,vertices);
    while true
        [minAng,worstFace] = minimumAngle(tri,(1:size(tri.ConnectivityList,1))');
        insertPt = tri.circumcenter(worstFace);
        [~,insertPt,f] = processPoint(tri.ConnectivityList,tri.Points,insertPt,tri.faceNormal, @distance_to_vertices,@distance_to_edges,@distance_to_surfaces, useSubSurface);        
        [~,r] = tri.incenter(f);
        [ tri2, is_new_face2, is_new_vertex2 ] = insert_vertex( tri, insertPt, f, r/2, is_new_face, is_new_vertex );
        [minAng2,~] = minimumAngle(tri2,(1:size(tri2.ConnectivityList,1))');
        if minAng2 > minAng+eps
            tri = tri2;
            is_new_face   = is_new_face2;
            is_new_vertex = is_new_vertex2;
        else
            break;
        end
    end

    faces2    = tri.ConnectivityList;
    vertices2 = tri.Points;
    new_faces_ID = find(is_new_face);
    
end

% return output arguments
distances      = D;  % (#qPoints x 1)
surface_points = P;  % (#qPoints x 3)
end


%% Non-vectorized Distance Functions
%  (can process only one point)


function [D,P,F] = processPoint(faces,vertices,point,normals, distance_to_vertices,distance_to_edges,distance_to_surfaces, useSubSurface)

d = NaN(3,1); % (distanceTypes x 1)
p = NaN(3,3); % (distanceTypes x xyz)
f = NaN(3,1); % (distanceTypes x 1)

% find nearest vertice
[d(1),p(1,:),f(1),v] = distance_to_vertices(faces,vertices,point,normals);
% d:  (1 x 1) signed distance to surface
% p:  (1 x 3) corresponding point on surface
% v:  (1 x 1) nearest vertex
% connectedFaces: (#connectedFaces x 1) face indices

if useSubSurface
    [tri2,~,faces_2To1] = subSurface( triangulation(faces,vertices), v, [], 2 );
    faces2 = tri2.ConnectivityList;
    vertices2 = tri2.Points;
    normals = normals(faces_2To1,:);
else
    faces2   = faces;
    vertices2 = vertices;
end

% find nearest point on all edges
[d(2),p(2,:),f(2)] = distance_to_edges(faces2,vertices2,point,normals);

% find nearest point on all surfaces
[d(3),p(3,:),f(3)] = distance_to_surfaces(faces2,vertices2,point,normals);

if useSubSurface
    % translate back f(2) and f(3)
    f(2:3) = faces_2To1(f(2:3));
end

% find minimum distance type
[~,I] = min(abs(d),[],1);
D = d(I);
P = p(I,:);
F = f(I);

end




function [D,P,F,V] = distance_to_vertices(faces,vertices,qPoint,normals)

% find nearest vertex
[D,nearestVertexID] = min(sum(bsxfun(@minus,vertices,qPoint).^2,2),[],1);
D = sqrt(D);
P = vertices(nearestVertexID,:); % (1 x 3)
V = nearestVertexID;

% find faces that belong to the vertex
connectedFaces = find(any(faces==nearestVertexID,2)); % (#connectedFaces x 1) face indices
assert(length(connectedFaces)>=1,'Vertex %u is not connected to any face.',nearestVertexID)
F = connectedFaces(1);
n = normals(connectedFaces,:); % (#connectedFaces x 3) normal vectors

% scalar product between distance vector and normal vectors
coefficients = dot2(n,qPoint-P);
sgn = signOfLargest(coefficients);
D = D*sgn;

end




function [D,P,F] = distance_to_edges(faces,vertices,qPoint,normals)

% Point-point representation of all edges
edges = [faces(:,[1,2]); faces(:,[1,3]); faces(:,[2,3])]; % (#edges x 2) vertice IDs

% Intersection between tangent of edge lines and query point
r1 = vertices(edges(:,1),:);   % (#edges x 3) first point of every edge
r2 = vertices(edges(:,2),:);   % (#edges x 3) second point of every edge
t = dot( bsxfun(@minus,qPoint,r1), r2-r1, 2) ./ sum((r2-r1).^2,2); % (#edges x 1) location of intersection relative to r1 and r2
t(t<=0) = NaN; % exclude intersections not between the two vertices r1 and r2
t(t>=1) = NaN;

% Distance between intersection and query point
P = r1 + bsxfun(@times,(r2-r1),t); % (#edges x 3) intersection points
D = bsxfun(@minus,qPoint,P); % (#edges x 3)
D = sqrt(sum(D.^2,2));       % (#edges x 1)
[D,I] = min(D,[],1);         % (1 x 1)
P = P(I,:);

% find faces that belong to the edge
inds = edges(I,:);  % (1 x 2)
inds = permute(inds,[3,1,2]);  % (1 x 1 x 2)
inds = bsxfun(@eq,faces,inds); % (#faces x 3 x 2)
inds = any(inds,3);    % (#faces x 3)
inds = sum(inds,2)==2; % (#faces x 1) logical indices which faces belong to the nearest edge of the query point
F = find(inds,1);
n = normals(inds,:); % (#connectedFaces x 3) normal vectors

% scalar product between distance vector and normal vectors
coefficients = dot2(n,qPoint-P); % (#connectedFaces x 1)
sgn = signOfLargest(coefficients);
D = D*sgn;

end




function [D,P,F] = distance_to_surfaces(faces,vertices,point,normals)

r1 = vertices(faces(:,1),:);   % (#faces x 3) % 1st vertex of every face
r2 = vertices(faces(:,2),:);   % (#faces x 3) % 2nd vertex of every face
r3 = vertices(faces(:,3),:);   % (#faces x 3) % 3rd vertex of every face

vq = bsxfun(@minus,point,r1);  % (#faces x 3)
D = dot(vq,normals,2);         % (#faces x 1) distance to surface
rD = bsxfun(@times,normals,D); % (#faces x 3) vector from surface to query point
P = bsxfun(@minus,point,rD);   % (#faces x 3) nearest point on surface; can be outside triangle

% find barycentric coordinates (query point as linear combination of two edges)
r31r31 = sum((r3-r1).^2,2);    % (#faces x 1)
r21r21 = sum((r2-r1).^2,2);    % (#faces x 1)
r21r31 = dot(r2-r1,r3-r1,2);   % (#faces x 1)
r31vq = dot(r3-r1,vq,2);       % (#faces x 1)
r21vq = dot(r2-r1,vq,2);       % (#faces x 1)

d = r31r31.*r21r21 - r21r31.^2;               % (#faces x 1)
bary = NaN(size(faces,1),3);                  % (#faces x 3)
bary(:,1) = (r21r21.*r31vq-r21r31.*r21vq)./d; % (#faces x 3)
bary(:,2) = (r31r31.*r21vq-r21r31.*r31vq)./d; % (#faces x 3)
bary(:,3) = 1 - bary(:,1) - bary(:,2);        % (#faces x 3)

% tri = triangulation(faces,vertices);
% bary = tri.cartesianToBarycentric((1:size(faces,1))',P); % (#faces x 3)

% exclude intersections that are outside the triangle
D( abs(d)<=eps | any(bary<=0,2) | any(bary>=1,2) ) = NaN;  % (#faces x 1)

% find nearest face for query point
[~,I] = min(abs(D),[],1); % (1 x 1)
D = D(I);       % (1 x 1)
P = P(I,:);     % (1 x 3)
F = I;          % (1 x 1)

end





function [tri2,highlight_faces2,faces_2To1,vertices_2To1] = subSurface( tri, vertex, highlight_faces, iterations )
% enlarge surface by vertex attachments

nFaces = size(tri.ConnectivityList,1);
nVertices = size(tri.Points,1);

faceL    = ind2logical([],nFaces);
verticeL = ind2logical([],nVertices);
verticeL(vertex) = true;
for x = 1:iterations
    faceIDs = tri.vertexAttachments(find(verticeL))'; %#ok<FNDSB> no array indexing 
    %faceIDs = cell2mat(faceIDs);
    faceIDs = [faceIDs{:}];
    faceL(faceIDs) = true;
    vertexIDs = tri.ConnectivityList(faceL,:);
    verticeL(vertexIDs) = true;
end

% faceL = any(tri.ConnectivityList==vertex,2);
% for x = 1:iterations-1
%     verts1 = tri.ConnectivityList(faceL,:);
%     faceL = any(ismember(tri.ConnectivityList,verts1),2);
% end

[connectedVertices,~,newIDs] = unique(tri.ConnectivityList(faceL,:));
faces2 = reshape(newIDs,[],3);               % (#connectedFaces    x 3) new face connectivity list
vertices2 = tri.Points(connectedVertices,:); % (#connectedVertices x 3) new vertice list

tri2 = triangulation(faces2,vertices2);

highlight_faces2 = ind2logical(highlight_faces,nFaces);
highlight_faces2 = highlight_faces2(faceL);
highlight_faces2 = find(highlight_faces2);

faces_2To1 = find(faceL);
vertices_2To1 = connectedVertices;

end


function logical = ind2logical(indices,len)
    logical = false(len,1);
    logical(indices(~isnan(indices))) = true;
end









%% Vectorized Distance Functions 
%  (can process more than one point on the same mesh)


function [D,P,F,V] = distance_to_vertices_vectorized(faces,vertices,points,normals)

% Requires Statistics Toolbox 
% [D,I] = pdist2(vertices,points, 'euclidean', 'Smallest',1); % (1 x #points)
% D = D'; % (#points x 1)
% I = I'; % (#points x 1)

D = sum(bsxfun(@minus,permute(vertices,[3,2,1]),points).^2,2); % (#points x 1 x #vertices)
[D,I] = min(D,[],3); % (#points x 1)
D = sqrt(D);         % (#points x 1)
P = vertices(I,:);   % (#points x 3)
V = I;

% find faces that belong to the vertex
inds = I; % (#points x 1)
inds = permute(inds,[3,2,1]); % (1 x 1 x #points)
inds = bsxfun(@eq,faces,inds); % (#faces x 3 x #points)
inds = any(inds,2); % (#faces x 1 x #points) logical indices which faces belong to the nearest edge of a query point
inds = permute(inds,[1,3,2]); % (#faces x #points)
inds = num2cell(inds,1); % (1 x #points) cell array with (#faces x 1) logical indices
n = cellfun(@(x) normals(x,:), inds, 'UniformOutput',false)'; % (#points x 1) cell array with (#connectedFaces x 3) normal vectors
F = cellfun(@(x) find(x,1), inds)'; % (#points x 1)

% scalar product between distance vector and normal vectors
coefficients = cellfun(@dot2, n, num2cell(points-P,2), 'UniformOutput',false);
sgn = cellfun(@signOfLargest,coefficients);
D = D.*sgn;

end



function [D,P,F] = distance_to_edges_vectorized(faces,vertices,points,normals)

euclid = @(A,dim) sqrt(sum(A.^2,dim));

% Point-point representation of all edges
edges = [faces(:,[1,2]); faces(:,[1,3]); faces(:,[2,3])]; % (#edges x 2) vertice IDs 

% Intersection between tangent of edge lines and query points
r1 = vertices(edges(:,1),:);   % (#edges x 3) first point of every edge 
r2 = vertices(edges(:,2),:);   % (#edges x 3) second point of every edge
qp = permute(points,[3,2,1]);  % (1 x 3 x #points) query points
t = bsxfun(@rdivide,...
        dot2(bsxfun(@minus,qp,r1), r2-r1),...
        sum((r2-r1).^2,2)); % (#edges x 1 x #points) location of intersection relative to r1 and r2 
t(t<=0) = NaN; % exclude intersections not between the two vertices r1 and r2  
t(t>=1) = NaN;

% Distance between intersection and query points
P = bsxfun(@plus,...
        r1,...
        bsxfun(@times,...
            (r2-r1),...
            t)); % (#edges x 3 x #points) intersection points
D = bsxfun(@minus,qp,P); % (#edges x 3 x #points) 
D = euclid(D,2);         % (#edges x 1 x #points) 
[D,I] = min(D,[],1);     % (1 x 1 x #points) 
D = squeeze(D);          % (#points x 1)
I = squeeze(I);          % (#points x 1)
P = permute(P,[2,1,3]);  % (3 x #edges x #points)
sz = [size(P) 1 1];
P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (3 x #points)
P = P';                  % (#points x 3)

% find faces that belong to the edge
inds = edges(I,:);  % (#points x 2)
inds = permute(inds,[4,3,1,2]); % (1 x 1 x #points x 2)
inds = bsxfun(@eq,faces,inds);  % (#faces x 3 x #points x 2)
inds = any(inds,4);    % (#faces x 3 x #points)
inds = sum(inds,2)==2; % (#faces x 1 x #points) logical indices which faces belong to the nearest edge of a query point
inds = permute(inds,[1,3,2]); % (#faces x #points)
inds = num2cell(inds,1);      % (1 x #points) cell array with (#faces x 1) logical indices
n = cellfun(@(x) normals(x,:), inds, 'UniformOutput',false)'; % (#points x 1) cell array with (#connectedFaces x 3) normal vectors
F = cellfun(@(x) find(x,1), inds)';  % (#points x 1)

% scalar product between distance vector and normal vectors
coefficients = cellfun(@dot2, n, num2cell(points-P,2), 'UniformOutput',false);
sgn = cellfun(@signOfLargest,coefficients);
D = D.*sgn;

end




function [D,P,F] = distance_to_surfaces_vectorized(faces,vertices,points,normals)

r1 = vertices(faces(:,1),:);   % (#faces x 3) % 1st vertex of every face 
r2 = vertices(faces(:,2),:);   % (#faces x 3) % 2nd vertex of every face 
r3 = vertices(faces(:,3),:);   % (#faces x 3) % 3rd vertex of every face 

qp = permute(points,[3,2,1]);  % (1 x 3 x #points) query points
vq = bsxfun(@minus,qp,r1);     % (#faces x 3 x #points) 
D = dot2(vq,normals);          % (#faces x 1 x #points) distance to surface
rD = bsxfun(@times,normals,D); % (#faces x 3 x #points) vector from surface to query point 
P = bsxfun(@minus,qp,rD);      % (#faces x 3 x #points) nearest point on surface; can be outside triangle 

% find barycentric coordinates (query point as linear combination of two edges) 
r31r31 = sum((r3-r1).^2,2);    % (#faces x 1)
r21r21 = sum((r2-r1).^2,2);    % (#faces x 1)
r21r31 = dot(r2-r1,r3-r1,2);   % (#faces x 1)
r31vq = dot2(r3-r1,vq);        % (#faces x 1 x #points)
r21vq = dot2(r2-r1,vq);        % (#faces x 1 x #points)

d = r31r31.*r21r21 - r21r31.^2; % (#faces x 1)
bary = NaN(size(faces,1), 2, size(points,1)); % (#faces x 3 x #points) 
bary(:,1,:) = bsxfun(@rdivide, bsxfun(@times,r21r21,r31vq) - bsxfun(@times,r21r31,r21vq), d); 
bary(:,2,:) = bsxfun(@rdivide, bsxfun(@times,r31r31,r21vq) - bsxfun(@times,r21r31,r31vq), d); 
bary(:,3,:) = 1 - bary(:,1,:) - bary(:,2,:);  % (#faces x 3 x #points) 

% exclude intersections that are outside the triangle
D( any(bary<=0,2) | any(bary>=1,2) ) = NaN;  % (#faces x 1 x #points)
D( abs(d)<=eps, :, : ) = NaN;

% find nearest face for every query point
[~,I] = min(abs(D),[],1); % (1 x 1 x #points)
I = squeeze(I); % (#points x 1)
D = D(sub2ind(size(D),I,ones(length(I),1),(1:length(I))'));
D = squeeze(D); % (#points x 1)
P = permute(P,[2,1,3]); % (3 x #faces x #points)
sz = [size(P) 1 1];
P = P(:,sub2ind(sz(2:3),I,(1:length(I))')); % (3 x #points)
P = P'; % (#points x 3)
F = I;  % (#points x 1)

end




function sgn = signOfLargest(coeff)
    [~,I] = max(abs(coeff));
    sgn = sign(coeff(I));
    if sgn==0, sgn=1; end
end




function d = dot2(A,B)
    % dot product along 2nd dimension with singleton extension
    d = sum(bsxfun(@times,A,B),2);
end









%% Insert Vertex into Triangulation
% 
% function [ tri2, is_new_face, is_new_vertex, new_vertex_ID ] = insert_vertex( tri, new_vertex, nearest_face, max_distance, is_new_face, is_new_vertex)
% 
% tri2 = tri;
% 
% % Adjacent triangles that have a dihedral angle that deviates from pi by an
% % angle greater than filterangle are preserved (no Delaunay flipping is performed) 
% % This conserves the previous triangulation shape. 
% filterangle = pi/3;
% 
% % Do nothing if new_vertex is close to any existing vertex
% nearest_face_def = tri2.ConnectivityList(nearest_face,:); % (1 x 3)   vertice IDs 
% nearest_face_vertices = tri2.Points(nearest_face_def,:);  % (3 x xyz) vertice coordinates  
% [distance,I] = min(sqrt(sum(bsxfun(@minus,new_vertex,nearest_face_vertices).^2,2)));
% nearest_vertex_ID = nearest_face_def(I);
% if distance<=max_distance
%     new_vertex_ID = nearest_vertex_ID;
%     return
% end
% 
% % Three cases:
% %   - Point is on a boundary edge
% %   - Point is on an edge inside the mesh
% %   - Point is inside a triangle (not on an edge)
% 
% % To create new faces, the old face is copied and vertice ID are replaced 
% % to maintain the normal orientation of the face. 
% 
% % Change to matrix representation to be able to edit the triangulation
% faces = tri2.ConnectivityList; % (#faces x 3)
% vertices = tri2.Points;        % (#vertices x 3)
% nVertices = size(tri2.Points,1);
% 
% % Distance to edge
% [edgeDist,pt] = distance_to_edges(nearest_face_def,vertices,new_vertex,tri2.faceNormal);
% 
% if ~isnan(edgeDist) && abs(edgeDist) <= max_distance 
%     % Point is on an edge
%     bary = tri2.cartesianToBarycentric(nearest_face,pt); % (1 x 3)
%     [~,I] = sort(bary); % ascending: first entry is the small one
%     faceDef = faces(nearest_face,:); % (1 x 3)
%     edge_verts = faceDef(I(2:3));    % (1 x 2) two vertices of edge
%     edge_faces = tri2.edgeAttachments(edge_verts); % (1 x 1) cell array 
%     edge_faces = edge_faces{1};      % (#edgeFaces x 1) faces at the edge
%     
%     FE = tri2.featureEdges(filterangle);
%     if ismember(edge_verts,FE,'rows') || ismember(flip(edge_verts),FE,'rows')
%         % Otherwise triangles with small angles will be produced
%         new_vertex = pt;
%     end
%         
%     if numel(edge_faces)==1 
%         % Boundary edge
%         
%         % Add 1 new vertex
%         vertices = [vertices;new_vertex];  % (#vertices x 3)
%         is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
%         new_vertex_ID = nVertices+1;       % scalar integer
%         
%         % Add 2 new faces
%         new_faces = [faces(nearest_face,:); faces(nearest_face,:)]; % (2 x 3) copy old face 
%         new_faces(1,I(2)) = new_vertex_ID; % (2 x 3) insert new vertex
%         new_faces(2,I(3)) = new_vertex_ID; % (2 x 3) insert new vertex
%         faces = [faces;new_faces];         % (#faces x 3)
%         is_new_face = [is_new_face;1;1];   % (#faces x 1)
%         
%         % Remove 1 old face
%         faces(nearest_face,:) = [];        % (#faces x 3)
%         is_new_face(nearest_face) = [];    % (#faces x 1)
%         
%         % Indices for Delaunay check
%         nFaces = size(faces,1);            % scalar integer
%         delaunay_check = nFaces-1:nFaces;  % (2 x 1) vector
%         
%     elseif numel(edge_faces)==2 
%     % Two attached faces on edge
%         
%         % Opposite vertices
%         edge_faces_def = faces(edge_faces,:);      % (2 x 3)
%         Lia = ismember(edge_faces_def,edge_verts); % logical (2 x 3), true where vertice IDs are on edge 
%         replace_face1 = find(Lia(1,:));            % (1 x 2) position of edge's vertice IDs in edge face 1  
%         replace_face2 = find(Lia(2,:));            % (1 x 2) position of edge's vertice IDs in edge face 2  
%         
%         % Add 1 new vertex
%         vertices = [vertices;new_vertex];  % (#vertices x 3)
%         new_vertex_ID = size(vertices,1);  % scalar integer
%         is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
%         
%         % Add 4 new faces
%         new_faces = [            % (4 x 3) copy old faces
%             edge_faces_def(1,:)
%             edge_faces_def(1,:)
%             edge_faces_def(2,:)
%             edge_faces_def(2,:)
%             ];
%         new_faces(1,replace_face1(1)) = new_vertex_ID; % (4 x 3) insert new vertex
%         new_faces(2,replace_face1(2)) = new_vertex_ID; 
%         new_faces(3,replace_face2(1)) = new_vertex_ID; 
%         new_faces(4,replace_face2(2)) = new_vertex_ID; 
%         
%         faces = [faces;new_faces];           % (#faces x 3)
%         is_new_face = [is_new_face;1;1;1;1]; % (#faces x 1)
%         
%         % Remove 2 old faces
%         faces(edge_faces,:) = [];     % (#faces x 3)
%         is_new_face(edge_faces) = []; % (#faces x 1)
%         
%         % Indices for Delaunay check
%         nFaces = size(faces,1);           % scalar integer
%         delaunay_check = nFaces-3:nFaces; % (1 x 4) vector
%         
%     else
%         error('More than two faces (%s) on an edge. This case is not implemented yet.',mat2str(edge_faces))
%     end    
% 
% else 
%     % Point is not on an edge
%     
%     % Add 1 new vertex
%     vertices = [vertices;new_vertex];  % (#vertices x 3)
%     is_new_vertex = [is_new_vertex;1]; % (#vertices x 1)
%     new_vertex_ID = nVertices+1;       % scalar integer
%     
%     % Add 3 new faces
%     nearest_face_def = faces(nearest_face,:); % (1 x 3)
%     new_faces = repmat(nearest_face_def,3,1); % (3 x 3) copy old face
%     new_faces(1,1) = new_vertex_ID;           % insert new vertex ID
%     new_faces(2,2) = new_vertex_ID;
%     new_faces(3,3) = new_vertex_ID;
%     faces = [faces;new_faces];         % (#faces x 3)
%     is_new_face = [is_new_face;1;1;1]; % (#faces x 1)
%     
%     % Remove 1 old face
%     faces(nearest_face,:) = [];     % (#faces x 1)
%     is_new_face(nearest_face) = []; % (#faces x 1)
%     
%     % Indices for Delaunay check
%     nFaces = size(faces,1);           % scalar
%     delaunay_check = nFaces-2:nFaces; % (1 x 3) vector
% end
% 
% % Change to triangulation class representation again
% tri2 = triangulation(faces,vertices);
% 
% % Restore Delaunay conditions
% for k = delaunay_check
%     [tri2, changed_faces] = makeDelaunay(tri2, k, filterangle);  
%     is_new_face = is_new_face | changed_faces; % (#faces x 1)
% end
% 
% end
% 
% 


%% Restore Delaunay Conditions
% 
% function [tri, changed_faces] = makeDelaunay( tri, face, filterangle )
% % Ensures that the Delaunay condition around a given face is met.
% % If an edge has to be flipped, the resulting changed faces are checked, too.
% % TODO Currently, more Delaunay checks than necessary are performed
% 
% changed_faces   = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces that have been flipped
% unchecked_faces = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces still have to be checked
% unchecked_faces(face) = true;
% 
% while any(unchecked_faces)
%     current_face = find(unchecked_faces,1);
%     unchecked_faces(current_face) = false;
%     [tri, changed_faces2] = makeDelaunay_neighbors( tri, current_face, filterangle );
%     unchecked_faces = unchecked_faces | changed_faces2;
%     changed_faces   = changed_faces   | changed_faces2;
% end
% 
% end
% 
% 
% 
% 
% function [tri, changed_faces] = makeDelaunay_neighbors( tri, face, filterangle )
% % Ensures that the Delaunay condition between the given face and its neighbors is met
% % After one flipping is performed, the corresponding faces are marked 
% % unchecked and the function returns. 
% 
% % For every neighbor face, check Delaunay condition and flip edge if necessary
% changed_faces = false(size(tri.ConnectivityList,1),1); % (#faces x 1) logical vector which faces that have been changed 
% 
% % Find neighbor faces and their vertices
% attached_faces = tri.neighbors(face)';      % (3 x 1) neighbor face IDs
% attached_faces(isnan(attached_faces)) = []; % (#attached_faces x 1) face IDs
% attached_vertice_IDs = tri.ConnectivityList(attached_faces,:); % (#attached_faces x 3) vertice IDs
% own_vertice_IDs = tri.ConnectivityList(face,:);       % (1 x 3) vertice IDs
% lia = ismember(attached_vertice_IDs,own_vertice_IDs); % (#attached_faces x 3) logical index
% attached_vertice_IDs_T = attached_vertice_IDs';       % (3 x #attached_faces) vertice IDs
% opposite_vertice_IDs = attached_vertice_IDs_T(~lia'); % (#attached_faces x 1) vertice IDs
% 
% for neighborID = 1:size(attached_faces,1)
%     face1 = face;
%     face2 = attached_faces(neighborID);
%     faces = [face1;face2]; % (2 x 1) face IDs
%     vertex1 = opposite_vertice_IDs(neighborID);
%     vertex2 = own_vertice_IDs(~ismember(own_vertice_IDs,attached_vertice_IDs(neighborID,:)));
%     vertices = [vertex1;vertex2]; % (2 x 1) vertice IDs
%     assert(length(vertices)==2)
%     if ~isDelaunay( tri, faces, vertices )
%         % Where sharp/acute/small angles occur, there is no Delaunay possible
%         % TODO if points to insert are outside the surface, the feature
%         % edges should be detected before any change of the mesh is performed
%         % (This is not the case here, since we use only surface projection points) 
%         FE = tri.featureEdges(filterangle);
%         edge = intersect(tri.ConnectivityList(faces(1),:),tri.ConnectivityList(faces(2),:));
%         sharpEdge = ismember(edge,FE,'rows') || ismember(flip(edge),FE,'rows');
%         % don't flip if other neighbors have more than one common vertex to flip partner 
%         neighbors = tri.neighbors(faces(1));
%         neighbors(isnan(neighbors)) = [];
%         neighbors(neighbors==faces(2)) = []; % neighbors of faces(1), without flip partner faces(2) 
%         nNeighbors = arrayfun(@(neighbor) length(intersect(tri.ConnectivityList(neighbor,:),tri.ConnectivityList(faces(2),:))), neighbors);
%         assert(~any(nNeighbors==0))
%         commonVertices = any(nNeighbors>1);
%         if ~sharpEdge && ~commonVertices
%             tri2 = flipFaces( tri, faces, vertices );
%             % Check if minimum angle increased by the flip
%             if minimumAngle(tri2,faces)>minimumAngle(tri,faces) 
%                 tri = tri2;
%                 changed_faces(faces) = true;
%                 return
%             end
%         end
%     end
% end
% end
% 
% 
% 
% 
% function [minAng,face, sortAngles,sortFaces] = minimumAngle( tri, faces )
% facesDef = tri.ConnectivityList(faces,:);   % (#faces x 3)
% facesCoord = tri.Points(facesDef(:),:);     % (#faces*3 x xyz)
% facesCoord = reshape(facesCoord,[],3,3);    % (#faces x 3 x xyz)
% ind = nchoosek(1:3,2);                      % (3 x 2)
% sideLength = facesCoord(:,ind(:,1),:)-facesCoord(:,ind(:,2),:);  % (#faces x 3 x xyz)
% sideLength = sqrt(sum(sideLength.^2,3));    % (#faces x 3)
% [~,r] = tri.circumcenter(faces);            % (#faces x 1)
% sinAngle = bsxfun(@rdivide,sideLength,2*r); % (#faces x 3)
% angle = asin(sinAngle);                     % (#faces x 3)
% minAngles = min(angle,[],2);                % (#faces x 1)
% [minAng,face] = min(minAngles);             % scalars
% [sortAngles,sortFaces] = sort(minAngles);   % (#faces x 1)
% assert(minAng==sortAngles(1))
% assert(face==sortFaces(1))
% end
% 
% 
% 
% 
% function isDelaunay = isDelaunay( tri, faces, vertices )
% % Checks if two faces (that share an edge) meet the Delaunay condition
% %    - tri:      triangulation object  
% %    - faces:    (2 x 1) face IDs
% %    - vertices: (2 x 1) vertice IDs of opposite vertices
% [CC,r] = tri.circumcenter(faces); % CC: (2 x xyz) centers;  r: (2 x 1) radius
% distance = sqrt(sum((tri.Points(vertices,:)-CC).^2,2)); % (2 x 1) distances
% isDelaunay = all( distance > r+eps );
% end
% 
% 
% 
% 
% function tri2 = flipFaces( tri, faces, vertices )
% tri2 = tri;
% % Change definition of two faces in the ConnectivityList
% face1_verts = tri2.ConnectivityList(faces(1),:);
% face2_verts = tri2.ConnectivityList(faces(2),:);
% chgVerts = face1_verts(face1_verts~=vertices(2));
% face1_verts(face1_verts==chgVerts(1)) = vertices(1);
% face2_verts(face2_verts==chgVerts(2)) = vertices(2);
% % Create updated triangulation object
% ConnectivityList = tri2.ConnectivityList;
% ConnectivityList(faces(1),:) = face1_verts;
% ConnectivityList(faces(2),:) = face2_verts;
% tri2 = triangulation(ConnectivityList,tri2.Points);
% end







