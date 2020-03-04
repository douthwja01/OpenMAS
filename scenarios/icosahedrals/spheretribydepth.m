function [vMat,fMat] = spheretribydepth(depth)
% SPHERETRIBYDEPTH is a high-performance vectorized function for building 
% a triangulation of a unit sphere based on recursive partitioning of each
% of Icosahedron faces into 4 triangles with vertices in the middles of 
% original face
%
% Input:
%   depth: double[1,1] - depth of partitioning, use 1 for the first level of
%       Icosahedron partitioning, and greater value for a greater level
%       of partitioning
%
% Output:
%   vMat: double[nVerts,3] - (x,y,z) coordinates of triangulation
%       vertices
%   fMat: double[nFaces,3] - indices of face verties in vertMat
%
% Example:
%   [vMat, fMat] = spheretribydepth(4);
%   patch('Vertices',vMat,'Faces',fMat,'FaceColor','g','EdgeColor','k');
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
if ~(numel(depth)&&isnumeric(depth)&&depth>=0&&fix(depth)==depth)
    error('spheretri:wrongInput',...
        'depth is expected to be a not negative integer scalar');
end
[vMat,fMat]=icosahedron();
[vMat,fMat]=shrinkfacetri(vMat,fMat,0,depth,@normvert);
end
function x=normvert(x)
x=x./repmat(realsqrt(sum((x.*x),2)),1,3);
end