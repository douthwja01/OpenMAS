function [vMat, fMat] = spheretri(nPoints)
% SPHERETRI is a high-performance vectorized function for building
% a triangulation of a unit sphere based on recursive partitioning of each
% of Icosahedron faces into 4 triangles with vertices in  the middles of 
% original face edgeMidMat. The function takes a number of  requested 
% points as an input. An exact number of points in the returned
% triangulation may not match this number, the function will choose a depth
% of Icosahedron partitioning minimally sufficient to provide the requested
% number of points
%
% Input:
%   nPoints: double[1,1] - number of requested points. An exact number of
%      points may not match this number, the function will choose a depth
%      of Icosahedron partitionin minimally sufficient for the requested
%      number of points
%
% Output:
%   vMat: double[nVerts,3] - (x,y,z) coordinates of triangulation
%       vertices
%   fMat: double[nFaces,3] - indices of face verties in vertMat
%
% Example:
%   [vMat, fMat] = spheretri(500);
%   patch('Vertices',vMat,'Faces',fMat,'FaceColor','g','EdgeColor','k');
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
if ~(isscalar(nPoints)&&isnumeric(nPoints)&&...
        (fix(nPoints)==nPoints)&&(nPoints>0))
    modgen.common.throwerror('wrongInput',...
        'nPoints is expected to be a positive integer scalar number');
end
partDepth=calcIcosPartDepth(nPoints);
[vMat, fMat] = spheretribydepth(partDepth);
%
function partDepth = calcIcosPartDepth( nPoints )
N_ICOSAHEDRON_VERTS=12;
N_ICOSAHEDRON_FACES=20;
N_ICOSAHEDRON_EDGES=30;
%
vertNum=N_ICOSAHEDRON_VERTS;
faceNum=N_ICOSAHEDRON_FACES;
edgeNum=N_ICOSAHEDRON_EDGES;
%
curDepth=0;
isStop=false;
while ~isStop
    curDepth=curDepth+1;
    vertNum=vertNum+edgeNum;
    edgeNum=2*edgeNum+3*faceNum;
    faceNum=4*faceNum;
    isStop=vertNum>=nPoints;
end
partDepth=curDepth;