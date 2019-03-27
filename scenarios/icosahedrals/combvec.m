function indMat = combvec(varargin)
% COMBVEC creates a matrix of combinations with elements from input vectors
%
% Usage: indMat=combvec(firstVec,secVec,...)
%
% Input:
%   optional:
%       firstVec: numeric[1,n1Elems] - first vector
%       ...
%       lastVec: numeric[1,nKElems] - last vector
%
% Output:
%   indMat: numeric[K,nCombs] - matrix with combinations from input vectors
%
% Example:
%
%     firstVec=[1,2,3];
%     secVec=[4,5];
%     indMat=mxberry.core.combvec(firstVec,secVec);
%     indMat
%
%     indMat =
%
%          1     2     3     1     2     3
%          4     4     4     5     5     5
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
nInputs=nargin;
if nInputs==0
    indMat = [];
else
    indMat = varargin{1};
    for iInp=2:nInputs
        curVec = varargin{iInp};
        indMat = [addBlockedInd(indMat,size(curVec,2));...
            addInterleavedInd(curVec,size(indMat,2))];
    end
end
%
function indResMat = addBlockedInd(indMat,nElems)
[nDims,nCombs] = size(indMat);
indResMat = zeros(nDims,nCombs*nElems);
indVec = 1:nCombs;
for iCombToAdd=(0:(nElems-1))*nCombs
    indResMat(:,indVec+iCombToAdd) = indMat;
end
%
function indResMat = addInterleavedInd(indMat,nElems)
[nDims,nCombs] = size(indMat);
indResMat = zeros(nDims*nElems,nCombs);
ind = 1:nDims;
for iCombToAdd=(0:(nElems-1))*nDims
    indResMat(ind+iCombToAdd,:) = indMat;
end
indResMat = reshape(indResMat,nDims,nElems*nCombs);