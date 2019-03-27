function [eMat,f2eMat,f2eIsDirMat] = mapface2edge(vMat,fMat)
% MAPFACE2EDGE creates a mapping from faces to edges 
%
% Input:
%   regular:
%       vMat: double[nVerts,3] - coordinates of vertices
%       fMat: double[nFace,3] - indices of face vertices in vMat  
%
% Output:  
%   eMat: double[nEdges,2] - contains indices of vertices corresponding
%       to each edge
%
%   f2eMat: double[nFaces,3] - contains indices of edges for
%       each face in this order (1-2, 2-3, 1-3)
%               
%   f2eIsDirMat: logical[nFaces,3] - contains true if face
%       references edge in a direct order (i.e. 1-2 for instance)
%       and false if reference is in an opposite order
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD, 
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
trObj = triangulation(fMat,vMat);
eMat=trObj.edges;
nEdges=size(eMat,1);
nFaces=size(fMat,1);
%% Build Face to Edges map and edge orientation map for each face
fMidCVec=trObj.edgeAttachments(eMat);
indShiftVec=cellfun('length',fMidCVec);
indF2EVec=zeros(sum(indShiftVec),1);
indF2EVec(1)=1;
indF2EVec(1+cumsum(indShiftVec(1:end-1)))=ones(nEdges-1,1);
indF2EVec=cumsum(indF2EVec);
indFVec=[fMidCVec{:}].';
%edges for each face are expected to be oriented as 
%1-2, 2-3 , 3-1
indEdgeNumVec=...
    all(fMat(indFVec,[1,2])==eMat(indF2EVec,:),2)...
    -all(fMat(indFVec,[2,1])==eMat(indF2EVec,:),2)...
    +2*all(fMat(indFVec,[2,3])==eMat(indF2EVec,:),2)...
    -2*all(fMat(indFVec,[3,2])==eMat(indF2EVec,:),2)...
    +3*all(fMat(indFVec,[1,3])==eMat(indF2EVec,:),2)...
    -3*all(fMat(indFVec,[3,1])==eMat(indF2EVec,:),2);
%
[~,indSortVec]=sortrows([indFVec,abs(indEdgeNumVec)]);
indF2EVec=indF2EVec(indSortVec);
f2eMat=reshape(indF2EVec,3,nFaces).';
f2eIsDirMat=reshape(indEdgeNumVec(indSortVec)>0,3,nFaces).';