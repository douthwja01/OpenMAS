function [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2, varargin)
%SURFACEINTERSECTION intersection of 2 surfaces
% [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2)
% calculates the intersection of surfaces 1 and 2. Code can either return
% just the matrix indicating which face of surface1 intersected with face
% of surface2, which is calculated using Tomas Moller algorithm, or can
% also return the actual line of intersection. In case when parts of the
% surface 1 and 2 lay on the same plane the intersection is a 2D area
% instead of 1D edge. In such a case the intersection area will be
% triangulated and intSurface.edges will hold the edges of the
% triangulation surface and intSurface.faces will hold the faces.
%
% INPUT:
%  * surface1 & surface2 - two surfaces defined as structs or classes.
%    Several inputs are possible:
%    - struct with "faces" and "vertices" fields
%    - 'triangulation' class (only the boundary surface will be used)
%    - 'delaunayTriangulation' class
%
% OUTPUT:
% * intMatrix - sparse Matrix with n1 x n2 dimension where n1 and n2 are
%               number of faces in surfaces
% * intSurface - a structure with following fields:
%     intSurface.vertices - N x 3 array of unique points
%     intSurface.edges    - N x 2 array of edge vertex ID's
%     intSurface.faces    - N x 3 array of face vertex ID's
%
% ALGORITHM:
% Based on Triangle/triangle intersection test routine by Tomas Möller, 1997.
%  See article "A Fast Triangle-Triangle Intersection Test",
%  Journal of Graphics Tools, 2(2), 1997
%  http://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf
%  http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt

%% Get FACES and VERTICES inputs
if isa(surface1, 'triangulation')
  [surface1.faces, surface1.vertices] = freeBoundary(surface1);
elseif isa(surface1, 'delaunayTriangulation')
  S = surface1;
  surface1 = [];
  surface1.faces    = S.ConnectivityList;
  surface1.vertices = S.Points;
  clear S
end
if isa(surface2, 'triangulation')
  [surface2.faces, surface1.vertices] = freeBoundary(surface2);
elseif isa(surface2, 'delaunayTriangulation')
  S = surface2;
  surface2 = [];
  surface2.faces    = S.ConnectivityList;
  surface2.vertices = S.Points;
  clear S
end
ok1 = isstruct(surface1) && isfield(surface1, 'vertices') && isfield(surface1, 'faces');
ok2 = isstruct(surface2) && isfield(surface2, 'vertices') && isfield(surface2, 'faces');
assert(ok1, 'Surface #1 must be a struct with "faces" and "vertices" fields' );
assert(ok2, 'Surface #2 must be a struct with "faces" and "vertices" fields' );

%% Flip dimentions if necessery
if size(surface1.faces,1)==3 && size(surface1.faces,2)~=3
  surface1.faces = surface1.faces';
end
if size(surface1.vertices,1)==3 && size(surface1.vertices,2)~=3
  surface1.vertices = surface1.vertices';
end
if size(surface2.faces,1)==3 && size(surface2.faces,2)~=3
  surface2.faces = surface2.faces';
end
if size(surface2.vertices,1)==3 && size(surface2.vertices,2)~=3
  surface2.vertices = surface2.vertices';
end

%% Parse extra parameters
getIntersection = (nargout>1);
debug = true;
PointRoundingTol = 1e6;
algorithm = 'moller';
k=1;
nVarargs = length(varargin);
while (k<=nVarargs)
  assert(ischar(varargin{k}), 'Incorrect input parameters')
  switch lower(varargin{k})
    case 'debug'
      debug = varargin{k+1}~=0;
      k = k+1;
    case 'algorithm'
      algorithm = lower(strtrim(varargin{k+1}));
      k = k+1;
    case 'pointroundingtol'
      PointRoundingTol = varargin{k+1};
      k = k+1;
  end
  k = k+1;
end

%% Initialize variables
epsilon = eps;
nFace1 = size(surface1.faces,1);
nFace2 = size(surface2.faces,1);
nVert1 = size(surface1.vertices,1);
nVert2 = size(surface2.vertices,1);

%% create strip down versions of MATLAB cross and dot function
cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));

%% Initialize output variables
% intersect is a nFace1 x nFace2 matrix. Possible values: -2 (do not know),
% -1 (coplanar with unknown overlap), 0 (no intersections), 1 (intersects).
% Negative values are internal only.
intMatrix  = zeros([nFace1,nFace2], 'int8')-2; % -2 indicates that there was no succesful test yet
intSurface.vertices = [];
intSurface.faces    = [];
intSurface.edges    = [];

% =======================================================================
%% === Stage 1 ==========================================================
% =======================================================================
% Each triangle is a subset of the plane it lies in, so for two triangles
% to intersect they must overlap along the line of intersection of their
% planes. Hence, a necessary condition for intersection is that each
% triangle must intersect the plane of the other.
% Möller’s method begins by checking the mutual intersection of each
% triangle with the plane of the other. To do so, it determines for each
% triangle on which side of the other triangle’s supporting plane its
% vertices lie. Now, if all vertices of one triangle lie on the same side
% and no vertex is on the plane, the intersection is rejected.

%% compute plane equations for each triangle of the surface #1
% plane equation #1: N1.X-d1=0
V1 = surface1.vertices(surface1.faces(:,1),:);
V2 = surface1.vertices(surface1.faces(:,2),:);
V3 = surface1.vertices(surface1.faces(:,3),:);
N1 = cross_prod(V2-V1,V3-V1); % array size nFace1 x 3
N1 = normalize(N1);
d1 = dot_prod(N1,V1);         % array size nFace1 x 1

%% Distance from surface #2 vertices to planes of surface #1
% Calculate signed distance from all vertices of surface #2 to each plane
% of of surface #1
du = zeros(nFace1,nVert2);
for iVert2 = 1:nVert2
  p = surface2.vertices(iVert2,:);
  du(:,iVert2) = N1(:,1)*p(1) + N1(:,2)*p(2) + N1(:,3)*p(3) - d1;
end
if debug
  assert(all(size(du)==[nFace1,nVert2]), 'Incorrect array dimensions: dv')
end
du(abs(du)<epsilon)=0; % robustness check
% Distances from vertex 1, 2 & 3 of faces of surface #2 to planes of surface #1
du1 = du(:,surface2.faces(:,1));
du2 = du(:,surface2.faces(:,2));
du3 = du(:,surface2.faces(:,3));
if debug
  assert(all(size(du1)==size(intMatrix)), 'Incorrect array dimensions: du1')
end
clear du
intMatrix(du1.*du2>0 & du1.*du3>0) = 0;   % same sign on all of them & not equal 0
if(all(intMatrix==0)), return; end        % no intersections
intMatrix(du1==0 & du2==0 & du3==0) = -1; % coplanar with unknown overlap

%% compute plane of triangle (U0,U1,U2)
% plane equation 2: N2.X-d2=0
U1 = surface2.vertices(surface2.faces(:,1),:);
U2 = surface2.vertices(surface2.faces(:,2),:);
U3 = surface2.vertices(surface2.faces(:,3),:);
N2 = cross_prod(U2-U1,U3-U1); % array size nFace1 x 3
N2 = normalize(N2);
d2 = dot_prod(N2,U1);        % array size nFace1 x 1

%% Distance from surface #1 vertices to planes of surface #2
% Calculate signed distance from all vertices of surface #1 to each plane
% of of surface #2
dv = zeros(nFace2,nVert1);
for iVert1 = 1:nVert1
  p = surface1.vertices(iVert1,:);
  dv(:,iVert1) = N2(:,1)*p(1) + N2(:,2)*p(2) + N2(:,3)*p(3) - d2;
end
if debug
  assert(all(size(dv)==[nFace2,nVert1]), 'Incorrect array dimensions: dv')
end
dv(abs(dv)<epsilon)=0; % robustness check
% Distances from vertex 1, 2 & 3 of faces of surface #1 to planes of surface #2
dv1 = dv(:,surface1.faces(:,1))';
dv2 = dv(:,surface1.faces(:,2))';
dv3 = dv(:,surface1.faces(:,3))';
if debug
  assert(all(size(dv1)==size(intMatrix)), 'Incorrect array dimensions: dv1')
end
clear dv
intMatrix(dv1.*dv2>0 & dv1.*dv3>0) = 0;   % same sign on all of them & not equal 0
if(all(intMatrix==0)), return; end        % no intersections
intMatrix(dv1==0 & dv2==0 & dv3==0) = -1; % coplanar with unknown overlap

% =======================================================================
%% === Stage 2 ==========================================================
% =======================================================================

%% Process remaining (non-coplanar) triangle pairs
tMsk = (intMatrix==-2);
n = nnz(tMsk);
if n>0
  [face1, face2] = find(tMsk);
  switch lower(algorithm)
    case 'moller'
      if size(dv1(tMsk),1)==1
        dv = [dv1(tMsk)', dv2(tMsk)', dv3(tMsk)'];
        du = [du1(tMsk)', du2(tMsk)', du3(tMsk)'];
      else
        dv = [dv1(tMsk), dv2(tMsk), dv3(tMsk)];
        du = [du1(tMsk), du2(tMsk), du3(tMsk)];
      end
      
      [intMatrix(tMsk), intSurface] = TriangleIntersection3D_Moller(...
        V1(face1,:), V2(face1,:), V3(face1,:), N1(face1,:), d1(face1,:), dv, ...
        U1(face2,:), U2(face2,:), U3(face2,:), N2(face2,:), d2(face2,:), du, ...
        getIntersection, debug);
    case 'rapid'
      % Undocumented experimental feature. In some experiments I got
      % identical results as with Moller algorithm, but others gave
      % different results. Often faster tham Moller.
      intMatrix(tMsk) = TriangleIntersection3D_Rapid( ...
        V1(face1,:), V2(face1,:), V3(face1,:), ...
        U1(face2,:), U2(face2,:), U3(face2,:), N1(face1,:), N2(face2,:) );
    otherwise
      error('Unknown algorithm name');
  end
end % if

%% Process coplanar triangle pairs. Pass #1:
% compare the overlap of the bounding boxes
tMsk = (intMatrix==-1);
if nnz(tMsk)>0
  [face1, face2] = find(tMsk);
  overlap = true;
  for idim = 1:3
    v = [V1(face1,idim), V2(face1,idim), V3(face1,idim)];
    u = [U1(face2,idim), U2(face2,idim), U3(face2,idim)];
    t1 = min(v,[],2);
    t2 = max(v,[],2);
    s1 = min(u,[],2);
    s2 = max(u,[],2);
    overlap = overlap & (s1<=t2 & t1<=s2);
  end
  % if overlap intMatrix will remain "-1" otherwise it will change to "0"
  intMatrix(tMsk) = -1*overlap;
  clear v u t1 t2 s1 s2 overlap
end

%% Process coplanar triangle pairs. Pass #2:
% use edge-edge intersections
tMsk = (intMatrix==-1);
if nnz(tMsk)>0
  [face1, face2] = find(tMsk);
  
  % repack data prior to function call
  V(:,:,1)=V1(face1,:); V(:,:,2)=V2(face1,:); V(:,:,3)=V3(face1,:);
  U(:,:,1)=U1(face2,:); U(:,:,2)=U2(face2,:); U(:,:,3)=U3(face2,:);
  [intMatrix(tMsk), intSurface2] = TriangleIntersection2D(V, U, ...
    N1(face1,:), getIntersection, debug);
  
  %% Merge surfaces
  if getIntersection
    np = size(intSurface.vertices,1);
    intSurface.vertices = [intSurface.vertices; intSurface2.vertices];
    intSurface.faces    = [intSurface.faces;    intSurface2.faces+np];
    intSurface.edges    = [intSurface.edges;    intSurface2.edges+np];
    if debug
      np = size(intSurface.vertices,1);
      assert(max(intSurface.faces(:))<=np, 'Bad surface definition')
      assert(max(intSurface.edges(:))<=np, 'Bad surface definition')
    end
  end
end

%% Clean up the outputs
intMatrix = sparse(double(intMatrix));
if(getIntersection)
  % make point array unique
  P = round(intSurface.vertices*PointRoundingTol)/PointRoundingTol;
  [~,ia,ic] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
  intSurface.vertices = intSurface.vertices(ia,:);
  intSurface.faces = ic(intSurface.faces);
  intSurface.edges = ic(intSurface.edges);
end
end % function

%% ========================================================================
function [iMsk, intSurface] = TriangleIntersection3D_Moller(...
  V1, V2, V3, N1, d1, dv, ...
  U1, U2, U3, N2, d2, du, ...
  getIntersection, debug)
%TriangleIntersection3D tests if 2 triangles defined in 3D intersect.
% This is a secondary test following Tomas Moller algorithm
%
% INPUTS:
%   V1, V2, V3, - Nx3 array of surface 1 triangle vertex coordinates
%   U1, U2, U3, - Nx3 array of surface 2 triangle vertex coordinates
%   N1, d1      - Nx3 array of surface 1 triangle plane equations N1.X-d1=0
%   N2, d2      - Nx3 array of surface 2 triangle plane equations N2.X-d2=0
%   dv          - Nx3 array of distances of surface 1 triangle vertices to surface 2 planes
%   du          - Nx3 array of distances of surface 2 triangle vertices to surface 1 planes
%   getIntersection - do we need to output the intersecting surface?
%      Algorithm is much simpler if we do not.
%   debug       - In the debugging mode much more extra "sanity check" test
%      are performed.
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%   intSurface - intersection surface
%
% ALGORITHM:
% The input triangles are guaranteed to intersect the line of intersection
% of the two planes. Furthermore, these intersections form intervals on
% this line, and the triangles overlap iff these intervals overlap as well.
% Hence, the last part of  the algorithm computes a parametric equation
% L(t) of the line of intersection of the two planes, finds the intervals
% (i.e. scalar intervals on L(t)) for which the line lies inside each
% triangle and performs a one-dimensional interval overlap test.
if debug
  ok = size(N1,2)==3 && size(N2,2)==3 && size(dv,2)==3 && size(du,2)==3 && ...
    size(V1,2)==3 && size(V2,2)==3 && size(V3,2)==3 && ...
    size(U1,2)==3 && size(U2,2)==3 && size(U3,2)==3;
  assert(ok, 'Incorrect array dimensions');
end

%% create strip down versions of MATLAB cross and dot function
cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
normalize = @(V) bsxfun(@rdivide,V, sqrt(sum(V.^2,2)));

%% Find intervals of surface 1 and 2 triangles
% compute the scalar intervals on L(t) for which the line lies inside each
% triangle

% Plane creates two open half-spaces. Find the odd vertex, which:
% 1) if no or two vertices are on the plane than pick the vertex which is
%    by itself in its half-space
% 2) if one vertex is on the plane and the other two occupy the same
%    half-space than pick the vertex on the plane
% 3) if one vertex is on the plane and the other two occupy different
%    half-spaces than pick one of the vertices off the plane
% Find vertex using a look-up table "lut" with key calculated based on
% sign of dv and du arrays
lut = [0;3;3;2;1;3;2;2;1;1;2;3;3;0;3;3;2;1;1;2;2;3;1;2;3;3;0];
n = numel(d1);
rows = (1:n)';

%% order surface 1 triangle vertices
a1 = lut(sign(dv)*[9; 3; 1] + 14); % calculate the key and call the look-up table
[b1, c1] = otherDim(a1);
if debug
  assert(all(a1>0), 'Something Wrong: triangles are coplanar')
end
a1 = sub2ind([n,3],rows,a1); % convert row and column IDs to array indecies
b1 = sub2ind([n,3],rows,b1);
c1 = sub2ind([n,3],rows,c1);

%% order surface 2 triangle vertices
a2 = lut(sign(du)*[9; 3; 1] + 14); % calculate the key and call the look-up table
[b2, c2] = otherDim(a2);
if debug
  assert(all(a2>0), 'Something Wrong: triangles are coplanar')
end
a2 = sub2ind([n,3],rows,a2);
b2 = sub2ind([n,3],rows,b2);
c2 = sub2ind([n,3],rows,c2);

%% compute direction of L the line of intersection of 2 planes
% containing 2 triangles. Line L parametric equation: t*D+O=0
D = cross_prod(N1,N2);    % D must be perpendicular to both N1 and N2
[~, maxDim] = max(abs(D),[],2); % compute and index to the largest component of D
if(getIntersection)
  D = normalize(D);
  O = zeros(n,3);
  d = [d1, d2, zeros(n,1)];
  for r =1:n
    N = [N1(r,:); N2(r,:); 0, 0, 0];
    N(3,maxDim(r)) = 1;
    dd = d(r,:)';
    O(r,:) = (N\dd)'; %Solve systems of linear equations N*D3 = d for D3
  end
  clear N d dd
end

%% projection of triangle(V1,V2,V3) and triangle(U1,U2,U3) onto intersection line
% Vp and Up are Nx3 arrays with columns indicating corners of triangles 1 and 2
if(getIntersection)
  Vp=[dot_prod(V1-O,D), dot_prod(V2-O,D), dot_prod(V3-O,D)];
  Up=[dot_prod(U1-O,D), dot_prod(U2-O,D), dot_prod(U3-O,D)];
else
  % Project on one of the axis (closest to the intersection line) instead.
  % Simplified projection is faster and sufficient if we do not need
  % intersection line
  idx = sub2ind([n,3],rows,maxDim);
  Vp = [V1(idx), V2(idx), V3(idx)];
  Up = [U1(idx), U2(idx), U3(idx)];
end
clear V1 V2 V3 U1 U2 U3

%% Calculate surface 1 and 2 triangle intervals
% t1 and t2 are intersection points of surface 1 with the intersection line
% t*D+O=0, and s1 & s2 are intersection points of surface 2 with the same
% line. Tomas Moller algorithm made this section much more complicated
% trying to avoid divisions. However, I could not detect any speed-up.
% Operations (ADD: 12; MUL:4 ; DIV:4 )
t1 = Vp(a1) - (Vp(b1)-Vp(a1)).*dv(a1)./(dv(b1)-dv(a1));
t2 = Vp(a1) - (Vp(c1)-Vp(a1)).*dv(a1)./(dv(c1)-dv(a1));
s1 = Up(a2) - (Up(b2)-Up(a2)).*du(a2)./(du(b2)-du(a2));
s2 = Up(a2) - (Up(c2)-Up(a2)).*du(a2)./(du(c2)-du(a2));

%% Order the intervals as to t1<t2 and s1<s2
msk = t2<t1; % order t1 and t2 so t1<t2
t = t1(msk); t1(msk)=t2(msk); t2(msk)=t; % swap
msk = s2<s1; % order s1 and s2 so s1<s2
t = s1(msk); s1(msk)=s2(msk); s2(msk)=t; % swap

%% Perform THE final test we were preparying for.
% It test for the overlap of 2 1D intervals s1->s2 and t1->t2
iMsk = (s1<t2 & t1<s2);

%% calculate intersection segments
n = nnz(iMsk);
if(getIntersection && n>0)
  % p1 = D*max(t1,s1) + O;    p2 = D*min(t2,s2) + O
  p1 = bsxfun(@times,D(iMsk,:),max(t1(iMsk),s1(iMsk))) + O(iMsk,:);
  p2 = bsxfun(@times,D(iMsk,:),min(t2(iMsk),s2(iMsk))) + O(iMsk,:);
  intSurface.vertices = [p1; p2];
  intSurface.faces    = [1:n; n+1:2*n; n+1:2*n]';
  intSurface.edges    = intSurface.faces(:,1:2);
else
  intSurface.vertices = [];
  intSurface.faces    = [];
  intSurface.edges    = [];
end % if
end % function

%% ========================================================================
function [overlap, intSurface] = TriangleIntersection2D(V, U, N, ...
  getIntersection, debug)
% Triangles V(V0,V1,V2) and U(U0,U1,U2) are are coplanar. Do they overlap?
% INPUTS:
% N - array(n,3) of surface normals where V(i,:,:) and U(i,:,:) are on the same plane
% V - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #1 vertices
% U - array(n,3,3) (nFace x 3 dimensions x 3 vertices) of surface #2 vertices
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%   intSurface - intersection surface


%  * parameters: vertices of triangle 1: V0,V1,V2
%  *             vertices of triangle 2: U0,U1,U2
%  * result    : returns 1 if the triangles intersect, otherwise 0

%% Constants needed for creating a mesh based on 3 to 6 points in a circle
tri_mesh{6}  = [1 2 6; 2 4 6; 2 3 4; 4 5 6];
tri_mesh{5}  = [1 2 3; 1 3 4; 4 5 1];
tri_mesh{4}  = [1 2 3; 1 3 4];
tri_mesh{3}  = 1:3;
vertices = [];
faces    = [];
pairs    = [];  % each row corresponds to pair of faces. match row number with face number
nVert    = 0;

%% use edge-edge intersections
overlap = false(size(N,1),1);
i1Idx = [1 1 1 2 2 2 3 3 3];
i2Idx = [3 3 3 1 1 1 2 2 2];
j1Idx = [1 2 3 1 2 3 1 2 3];
j2Idx = [3 1 2 3 1 2 3 1 2];
for row = 1:size(N,1)
  % When it is necesary to project 3D plane on 2D, dIdx will be the optimal
  % dimensions to use.
  [~, a] = max(abs(N(row,:))); 
  [b, c] = otherDim(a); 
  dIdx = [b, c]; 
  order = [];

  %% test all edges of triangle 1 against the edges of triangle 2
  % triangles overlap if edges cross
  [edgeMat, P] = EdgesIntersect3D(...
    squeeze(V(row,:,i1Idx))',squeeze(V(row,:,i2Idx))', ...
    squeeze(U(row,:,j1Idx))',squeeze(U(row,:,j2Idx))');
  overlap(row) = any(edgeMat);
  if ~getIntersection && overlap(row), continue; end
  
  if ~overlap(row)
    %% project onto an axis-aligned plane, that maximizes the area
    % of the triangles, compute indices: dIdx which correspond to 2 smallest N1
    % components.
    V2d = [V(row,dIdx,1); V(row,dIdx,2); V(row,dIdx,3)]; % each row is a 2D vertex
    U2d = [U(row,dIdx,1); U(row,dIdx,2); U(row,dIdx,3)];
    
    %% test if tri1 is totally contained in tri2 or vice varsa
    if PointInTriangle2D(V2d(1,:), U2d) % tri1 is totally contained in tri2
      overlap(row) = true;
      order = 1:3;
    elseif PointInTriangle2D(U2d(1,:), V2d) % tri2 is totally contained in tri1
      overlap(row) = true;
      order = 4:6;
    end
    if overlap(row) && ~getIntersection, continue; end
    clear V2d U2d
  end
  
  %% Build the intersection surface
  if getIntersection && overlap(row)
    %Assemble all the points which might be needed for desining
    %intersection polygon: Intersection points and points from triangle 1
    %and 2
    points   = [P(edgeMat,:); squeeze(V(row,:,:))'; squeeze(U(row,:,:))'];
    if isempty(order) % when one tri is totally contained in the other tri then order is set
      order = IntersectionPolygon(edgeMat>0, points, dIdx, debug);
      if isempty(order), continue; end
    end
    nPoint   = length(order);    % how many points will be added?
    nFace    = nPoint-2;         % how many faces will be added?
    vertices = [vertices; points(order,:)]; %#ok<*AGROW>
    faces    = [faces; nVert+tri_mesh{nPoint} ];
    pairs    = [pairs; row+zeros(nFace,1)];  % each row corresponds to pair of faces. match row number with face number
    nVert    = nVert + nPoint;
    if debug
      assert(max(faces(:))<=size(vertices,1), 'Bad surface definition')
    end
  end
end % for

%% Prepare outputs
intSurface.vertices = vertices;
intSurface.faces    = faces;
if isempty(faces)
  intSurface.edges = [];
else
  intSurface.edges = [faces(:,1:2); faces(:,2:3); faces(:,[1,3])];
end
end % function

%% ========================================================================
function polygon = IntersectionPolygon(edgeMat, points, dIdx, debug)
% edgeMat is an edge intersection matrix with 3 rows for edges between
% the points 1-3, 1-2, & 2-3 of the triangle 1 and 3 columns for the same
% edges of the triangle 2. If 2 edges intersect a point of intersection
% is calculated and stored in array "points" followed by points of the
% triangles 1 & 2.  This function calculates the polygon of the intersection
% between 2 triangles.

persistent orderLUT verified
if isempty(orderLUT) || isempty(orderLUT{3})
  % This pre-calculated look-up table is used to quickly look up the order of
  % the vertices in array "points" which make up polygon of the intersection
  % between 2 triangles. A unique key is calculated for each edgeMat using
  % dot product between edgeMat(:) and [256 128 64 32 16 8 4 2 1], which is
  % used to look up point order around the polygon. Negative numbers in the
  % LUT indicate values which were not observed yet so they were not
  % independently verified.
  % reshape(sprintf('%09s',dec2base(key, 2)),3,3) will convert from the key
  % to matrix.
  OrderLUT = zeros(432,1);  
  OrderLUT(003) = 127;
  OrderLUT(005) = 128;
  OrderLUT(006) = 126;
  OrderLUT(009) = 124;
  OrderLUT(010) = 1427;
  OrderLUT(012) = 1428;
  OrderLUT(017) = 1427;
  OrderLUT(018) = 124;
  OrderLUT(020) = 1426;
  OrderLUT(024) = 127;
  OrderLUT(027) = 1243;
  OrderLUT(029) = 12438;
  OrderLUT(030) = 12034;
  OrderLUT(033) = 1428;
  OrderLUT(034) = 1426;
  OrderLUT(036) = 124;
  OrderLUT(040) = 128;
  OrderLUT(043) = 21834;
  OrderLUT(045) = 1243;
  OrderLUT(046) = 21349;
  OrderLUT(048) = 126;
  OrderLUT(051) = 12340;
  OrderLUT(053) = 12943;
  OrderLUT(054) = 1243;
  OrderLUT(065) = 125;
  OrderLUT(066) = 1527;
  OrderLUT(068) = 1825;
  OrderLUT(072) = 123;
  OrderLUT(080) = 1327;
  OrderLUT(083) = 15234;
  OrderLUT(085) = -15234;
  OrderLUT(086) = -15243;
  OrderLUT(090) = 13247;
  OrderLUT(092) = -13247;
  OrderLUT(096) = 1328;
  OrderLUT(099) = 152834;
  OrderLUT(101) = 15234;
  OrderLUT(102) = 152349;
  OrderLUT(106) = 132847;
  OrderLUT(108) = 13247;
  OrderLUT(114) = 102347;
  OrderLUT(116) = -13247;
  OrderLUT(129) = 1527;
  OrderLUT(130) = 125;
  OrderLUT(132) = 1526;
  OrderLUT(136) = 1327;
  OrderLUT(139) = 15243;
  OrderLUT(141) = 152438;
  OrderLUT(142) = 152034;
  OrderLUT(144) = 123;
  OrderLUT(153) = 12347;
  OrderLUT(156) = 123047;
  OrderLUT(160) = 1326;
  OrderLUT(163) = -152043;
  OrderLUT(165) = 13247;
  OrderLUT(166) = 15234;
  OrderLUT(169) = -182347;
  OrderLUT(172) = 193247;
  OrderLUT(177) = -132047;
  OrderLUT(180) = 13247;
  OrderLUT(192) = 127;
  OrderLUT(195) = 1243;
  OrderLUT(197) = 12438;
  OrderLUT(198) = 12034;
  OrderLUT(202) = 12364;
  OrderLUT(204) = 123648;
  OrderLUT(209) = 21364;
  OrderLUT(212) = -21364;
  OrderLUT(216) = 1243;
  OrderLUT(225) = -124638;
  OrderLUT(226) = 120364;
  OrderLUT(232) = 12438;
  OrderLUT(238) = 124356;
  OrderLUT(240) = 12034;
  OrderLUT(245) = -214356;
  OrderLUT(257) = 1528;
  OrderLUT(258) = 1526;
  OrderLUT(260) = 125;
  OrderLUT(264) = 1328;
  OrderLUT(267) = -152438;
  OrderLUT(269) = 15243;
  OrderLUT(270) = -152943;
  OrderLUT(272) = 1326;
  OrderLUT(275) = 152340;
  OrderLUT(277) = 152943;
  OrderLUT(278) = 15243;
  OrderLUT(281) = 182347;
  OrderLUT(282) = -103247;
  OrderLUT(288) = 123;
  OrderLUT(297) = 12347;
  OrderLUT(298) = -123947;
  OrderLUT(305) = 123947;
  OrderLUT(306) = 12347;
  OrderLUT(320) = 128;
  OrderLUT(323) = 21834;
  OrderLUT(325) = 1243;
  OrderLUT(326) = 21349;
  OrderLUT(330) = -123648;
  OrderLUT(332) = 12364;
  OrderLUT(337) = 183642;
  OrderLUT(340) = -129364;
  OrderLUT(344) = 21834;
  OrderLUT(350) = -124365;
  OrderLUT(353) = 12463;
  OrderLUT(354) = 136492;
  OrderLUT(360) = 1243;
  OrderLUT(368) = 12943;
  OrderLUT(371) = 126543;
  OrderLUT(384) = 126;
  OrderLUT(387) = 12340;
  OrderLUT(389) = 12943;
  OrderLUT(390) = 1243;
  OrderLUT(394) = -103642;
  OrderLUT(396) = 129364;
  OrderLUT(401) = 123640;
  OrderLUT(404) = 12364;
  OrderLUT(408) = 12340;
  OrderLUT(413) = 215643;
  OrderLUT(417) = -136492;
  OrderLUT(418) = 12463;
  OrderLUT(424) = 13492;
  OrderLUT(427) = -213456;
  OrderLUT(432) = 1342;
  
  % Convert to more convinient format
  orderLUT = cell(size(OrderLUT));
  for i = 1:size(OrderLUT,1)
    polygon = abs(OrderLUT(i));
    if polygon>0
      polygon = num2str(polygon)-48; % Convert from a single number to array of digits
      polygon(polygon==0) = 10;      % 0 stands for 10
      orderLUT{i} = polygon;
    end
  end
  % Negative numbers in the LUT indicate values which were not observed yet
  % so they were not independently verified.
  verified = OrderLUT>0;
  clear OrderLUT
end

%% Calculate unique key for each edgeMat configuration
key = dot(1*edgeMat(:)', [256 128 64 32 16 8 4 2 1]);
assert(key<=432, 'Error: in IntersectionPolygon: key is out of bound');

%% Look up the point order around the polygon
polygon = orderLUT{key};
if (isempty(polygon))
  return
end

%% in a rare case of 2 intersections there is ambiguity if one or two
% vertices of the triangle lay inside the other triangle. OrderLUT stores
% only the single vertex cases.
nx = nnz(edgeMat(:));
if nx==2
  pList = polygon;       % list of vertices to check
  pList(pList<=nx) = []; % keep only the triangle points of the polygon
  flip = false;    % was there a flip from single vertex to vertices case?
  for ip = 1:length(pList)
    p = pList(ip);                 % point to check
    t = floor((p-nx-1)/3);         % does it belong to triangle 0 or 1 (actually 1 or 2)
    tri = (1:3) + nx + 3*abs(1-t); % Points belonging to the other triangle
    if ~PointInTriangle2D(points(p,dIdx), points(tri,dIdx))
      d = nx+t*3;    % offset
      % "p-d" is vertex number of point just tested: 1, 2, or 3. "b, c" are
      % the other 2 vertices
      [b, c] = otherDim(p-d);
      polygon = [polygon(polygon~=p), b+d, c+d]; % remove i2 and add i0 and i1
      flip = true;
    end
  end
  if flip
    % if ther were any flips than use existing codes to figure out the
    % order of the points around the polygon
    DT = delaunayTriangulation(points(polygon,dIdx));
    idx = freeBoundary(DT)';
    idx(2,:) = [];
    polygon = polygon(idx);
  end
end

%% Check to duplicate points
tol = 1e6;
P = round(points(polygon,:)*tol)/tol;
[~,ia] = unique(P,'rows'); % V = P(ia,:) and P = V(ic,:).
polygon = polygon(sort(ia));

%% Test the results using more expensive function
doPlot = (~verified(key));
if debug && length(polygon)>3
  DT = delaunayTriangulation(points(polygon,dIdx));
  idx = freeBoundary(DT)';
  idx(2,:) = [];
  k = max(abs(diff(idx)));
  %doPlot = (k>1 && k<(length(idx)-1)) || (~verified(key));
  assert(k==1 || k==(length(idx)-1), 'Two triangle intersection polygon is not convex')
end
if debug && doPlot % plot the interesting cases
  PlotTwoTriangles(points, polygon, 'm')
  title(sprintf('key = %i', key));
end 

end % function

%% ========================================================================
function PlotTwoTriangles(points, polygon, color)
% Plotting function used for debugging
nx = size(points,1)-6;
d = (max(points,[],1)-min(points,[],1))/200;
figure(2)
clf
hold on
line( points(nx+(1:2),1), points(nx+(1:2),2), points(nx+(1:2),3), 'Color', 'g');
line( points(nx+(2:3),1), points(nx+(2:3),2), points(nx+(2:3),3), 'Color', 'g');
line( points(nx+[1,3],1), points(nx+[1,3],2), points(nx+[1,3],3), 'Color', 'g');
line( points(nx+(4:5),1), points(nx+(4:5),2), points(nx+(4:5),3), 'Color', 'b');
line( points(nx+(5:6),1), points(nx+(5:6),2), points(nx+(5:6),3), 'Color', 'b');
line( points(nx+[4,6],1), points(nx+[4,6],2), points(nx+[4,6],3), 'Color', 'b');
plot3( points(:,1), points(:,2), points(:,3), 'm.');
if (length(polygon)>2)
  idx = polygon([1:end, 1]);
  plot3( points(idx,1), points(idx,2),points(idx,3), 'Color', color, 'LineWidth', 1);
end
for i = 1:nx+6
  text(points(i,1)+d(1), points(i,2)+d(2), points(i,3), num2str(i))
end

end % function

%% ========================================================================
function [intersect, X] = EdgesIntersect3D(V1,V2, U1,U2)
%EdgesIntersectPoint3D calculates point of intersection of 2 coplanar
% segments in 3D
%
% INPUTS:
%   V1,V2 - 1 x 3 coordinates of endpoints of edge 1
%   U1,U2 - 1 x 3 coordinates of endpoints of edge 2
% OUTPUT:
%   X - 1 x 3 coordinates of the intersection point
A = V2-V1;
B = U1-U2;
C = U1-V1;
%% Solve system of equations [A,B,1] * [d;e;0] = C for d and e
det3 = @(a,b) ... % determinant of a matrix with columns: [a, b, 1]
  a(:,1).*b(:,2)-a(:,3).*b(:,2) + ...
  a(:,2).*b(:,3)-a(:,2).*b(:,1) + ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3);
f=det3(A,B); % https://en.wikipedia.org/wiki/Cramer%27s_rule#Explicit_formulas_for_small_systems
t=det3(C,B)./f; % use Cramer's rule
s=det3(A,C)./f;
intersect = (t>=0 & t<=1 & s>=0 & s<=1);
X = V1 + bsxfun(@times,A,t);
end % function

%% ========================================================================
function inside = PointInTriangle2D(V1, U)
% check if V1 is inside triangle U (U1,U2,U3)
% Algorithm is checking on which side of the half-plane created by the
% edges the point is. It uses sign of determinant to calculate orientation
% of point triplets.
% INPUTS:
%   V1 - 1 x 2 coordinates of a point
%   U  - 3 x 2 coordinates of endpoints of 3 edges of a triangle
% OUTPUT:
%   inside - a boolean or boolean array
det2 = @(A,B,C) (A(:,1)-C(:,1))*(B(:,2)-C(:,2)) - (B(:,1)-C(:,1))*(A(:,2)-C(:,2));
b1 = (det2(U(1,:), U(2,:), V1) > 0);
b2 = (det2(U(2,:), U(3,:), V1) > 0);
b3 = (det2(U(3,:), U(1,:), V1) > 0);
inside = ((b1 == b2) & (b2 == b3)); % inside if same orientation for all 3 edges
end % function

%% ========================================================================
function [b, c] = otherDim(a)
% return set [1 2 3] without k
b = mod(a+1,3)+1;  % b and c are vertices which are on the same side of the plane
c = 6-a-b;         % a+b+c = 6
end


%% ========================================================================
function overlap = TriangleIntersection3D_Rapid( v1, v2, v3, u1, u2, u3, n1, n2 )
%TriangleIntersection3D tests if 2 triangles defined in 3D intersect.
%
% INPUTS:
%   v1, v2, v3, - Nx3 array of surface 1 triangle vertex coordinates
%   u1, u2, u3, - Nx3 array of surface 2 triangle vertex coordinates
%   n1, n2      - Nx3 array of surface 1 & 2 triangle plane normals. Those
%      are optional and if provided than the first 2 steps of the algorithm
%      (which are equivalent to first 2 steps of Moller algorithm) will be 
%      skipped.
%
% OUTPUT:
%   iMsk - N x 1 intersection boolean mask marking which triangles overlap
%
% ALGORITHM:
%   translated from the UNC-CH V-Collide RAPID code
%    https://wwwx.cs.unc.edu/~geom/papers/COLLISION/vcol.pdf

global V1 V2 V3 U1 U2 U3

cross_prod = @(a,b) [...
  a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  a(:,1).*b(:,2)-a(:,2).*b(:,1)];

%% shift t1 and t2 by p1#
V1 = zeros(size(v1));
V2 = v2-v1;
V3 = v3-v1;
U1 = u1-v1;
U2 = u2-v1;
U3 = u3-v1;
clear v1 v2 v3 u1 u2 u3

if(nargin<7)
  %% now begin the series of tests
  n1 = cross_prod( V2-V1, V3-V1 ); % face normals
  n2 = cross_prod( U2-U1, U3-U1 ); % face normals
end
  
%% test the face normals
overlap = project6(n1) & project6(n2);
V1 = V1(overlap,:);
V2 = V2(overlap,:);
V3 = V3(overlap,:);
U1 = U1(overlap,:);
U2 = U2(overlap,:);
U3 = U3(overlap,:);
n1 = n1(overlap,:);
n2 = n2(overlap,:);

%% compute triangle edges
e1 = V2-V1;
e2 = V3-V2;
e3 = V1-V3;
f1 = U2-U1;
f2 = U3-U2;
f3 = U1-U3;

%% run more tests
overlap2 = project6(cross_prod(e1, f1));
overlap2 = project6(cross_prod(e1, f2)) & overlap2;
overlap2 = project6(cross_prod(e1, f3)) & overlap2;
overlap2 = project6(cross_prod(e2, f1)) & overlap2;
overlap2 = project6(cross_prod(e2, f2)) & overlap2;
overlap2 = project6(cross_prod(e2, f3)) & overlap2;
overlap2 = project6(cross_prod(e3, f1)) & overlap2;
overlap2 = project6(cross_prod(e3, f2)) & overlap2;
overlap2 = project6(cross_prod(e3, f3)) & overlap2;
overlap2 = project6(cross_prod(e1, n1)) & overlap2;
overlap2 = project6(cross_prod(e2, n1)) & overlap2;
overlap2 = project6(cross_prod(e3, n1)) & overlap2;
overlap2 = project6(cross_prod(f1, n2)) & overlap2;
overlap2 = project6(cross_prod(f2, n2)) & overlap2;
overlap2 = project6(cross_prod(f3, n2)) & overlap2;
overlap(overlap) = overlap2;
end

%% ========================================================================
function pass = project6( p )
% project all 6 vertices of both triangles onto vector p and check if two
% projections overlap
global V1 V2 V3 U1 U2 U3
dot_prod = @(a,b) a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
%% Project vertices of triangle 1 and find the bounds min1 and max1
P = [dot_prod(p, V1), dot_prod(p, V2), dot_prod(p, V3)];
max1 = max(P,[],2);
min1 = min(P,[],2);
%% Project vertices of triangle 2 and find the bounds min1 and max1
P = [dot_prod(p, U1), dot_prod(p, U2), dot_prod(p, U3)];
max2 = max(P,[],2);
min2 = min(P,[],2);
%% Compare the bounds to see if they overlap
pass = (( min1 < max2 ) & ( min2 < max1 )) | ~dot_prod(p, p);
end
