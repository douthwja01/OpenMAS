function [vMat,fMat,SStats,eMat,f2eMat,f2eIsDirMat]=shrinkfacetri(vMat,...
    fMat,maxEdgeLength,nMaxSteps,fVertAdjustFunc)
% SHRINKFACETRI shrinks faces of 3D triangulation space down to a
% degree where a length of each edge is less or equal to a specified value
%
% Input:
%   regular:
%       vMat: double[nVerts,3] - coordinates of vertices
%       fMat: double[nFace,3] - indices of face vertices in vMat
%       maxEdgeLength: double[1,1] - maximum allowed edge length in the
%          resulting triangulation
%
%   optional:
%       nMaxSteps: double[1,1] - maximum allowed step number - edge
%          shrinking is an iterative process so this parameter limits
%          a number of steps, by default the number of steps is unlimited
%          (nMaxSteps=Inf)
%       fVertAdjustFunc: function_handle[1,1] - function responsible for
%          transforming vertices on each step, the function accept vMat and
%          output vMat, by default the function is @deal (i.e. no
%          transformation is performed)
% Output:
%       vMat: double[nNewVerts,3] - coordinates of vertices in the
%           resulting triangulation
%       fMat: double[nNewFace,3] - indices of face vertices in vMat
%           indices of face vertices in the resulting triangulation
%
%       SStats: struct[1,1] - structure containing the statistics
%           related to the face shrink process, includes the following
%           fields:
%               nSteps: double[1,1] - number of performed steps
%               nVertVec: double[nSteps+1,1] - number of vertices on each
%                  step
%               nFaceVec: double[nSteps+1,1] - numbers of faces on each 
%                   step, the first item corresponds to zero step
%               nEdgeVec: double[nSteps+1,1] - numbers of edges on each 
%                   step
%               nEdgesToShrinkVec: double[nSteps+1,1] - numbers of edges to
%                  shrink on each step
%               maxEdgeLengthVec: double[nSteps+1,1] - vector of maximum
%                  edge lengths on each step
%
%       eMat: double[nEdges,2] - contains indices of vertices corresponding
%          to each edge
%
%       f2eMat: double[nFaces,3] - contains indices of edges for
%           each face in this order (1-2, 2-3, 1-3)
%               
%       f2eIsDirMat: logical[nFaces,3] - contains true if face
%           references edge in a direct order (i.e. 1-2 for instance)
%           and false if reference is in an opposite order
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD, 
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
if nargin<5
    isAdjustFuncSpec=false;
    if nargin<4
        nMaxSteps=Inf;
    end
else
    isAdjustFuncSpec=true;
end
%
isStatCollected=nargout>2;
%
%% Build Face to Edges map and edge orientation map for each face
[eMat,f2eMat,f2eIsDirMat]=mapface2edge(vMat,fMat);
nVerts=size(vMat,1);
nEdges=size(eMat,1);
nFaces=size(fMat,1);
%
%% Build edge distances
dMat=vMat(eMat(:,1),:)-vMat(eMat(:,2),:);
eLengthVec=realsqrt(sum(dMat.*dMat,2));
%% Set up statistics collection
if isStatCollected
    nVertVec=nVerts;
    nEdgeVec=nEdges;
    nFaceVec=nFaces;
    nEdgesToShrinkVec=zeros(0,1);
    maxEdgeLengthVec=zeros(0,1);
end
%%
iStep=1;
while true
    %% Find edges that needs to be shortened
    isE2OrigPartVec=eLengthVec>maxEdgeLength;
    if any(isE2OrigPartVec)&&(iStep<=nMaxSteps)
        %% Find faces that needs to be shortened
        if nFaces>1
            isF2PartVec=any(isE2OrigPartVec(f2eMat),2);
        else
            %this is because isE2PartVec(f2eMat) produces a column-vector
            %for one face
            isF2PartVec=any(isE2OrigPartVec(f2eMat));
        end
        %% Readjust indices of partitioned edges
        f2ePartMat=f2eMat(isF2PartVec,:);
        f2eIsDirPartMat=f2eIsDirMat(isF2PartVec,:);
        isE2PartVec=isE2OrigPartVec;
        isE2PartVec(f2ePartMat)=true;
        %% Build indices of partitioned vertices
        indV1Vec=eMat(isE2PartVec,1);
        indV2Vec=eMat(isE2PartVec,2);
        vNewMat=(vMat(indV1Vec,:)+vMat(indV2Vec,:))*0.5;
        nNewVerts=size(vNewMat,1);
        %% Collect stats
        if isStatCollected
            nEdgesToShrinkVec=[nEdgesToShrinkVec;sum(isE2PartVec)]; %#ok<*AGROW>
            maxEdgeLengthVec=[maxEdgeLengthVec;max(eLengthVec)];
        end
        %% Find new faces, edges and vertices
        nShrinkedFaces=sum(isF2PartVec);
        %% Build indices of new vertices for each partitioned face
        indEPartVec=cumsum(isE2PartVec)+nVerts;
        indVF12Vec=indEPartVec(f2ePartMat(:,1));
        indVF23Vec=indEPartVec(f2ePartMat(:,2));
        indVF13Vec=indEPartVec(f2ePartMat(:,3));
        %
        %% Build indices of vertices of partitiened faces
        %
        indVF1Vec=fMat(isF2PartVec,1);
        indVF2Vec=fMat(isF2PartVec,2);
        indVF3Vec=fMat(isF2PartVec,3);
        %% Remove partitioned faces and edges
        fMat=fMat(~isF2PartVec,:);
        f2eMat=f2eMat(~isF2PartVec,:);
        f2eIsDirMat=f2eIsDirMat(~isF2PartVec,:);
        isEKeptVec=false(nEdges,1);
        isEKeptVec(f2eMat)=true;
        indEKeptVec=cumsum(isEKeptVec);
        if size(f2eMat,1)==1
            %this is because f2eMat becomes a column for 1 face in f2eMat
            f2eMat=transpose(indEKeptVec(f2eMat));
        else
            f2eMat=indEKeptVec(f2eMat);
        end
        eMat=eMat(isEKeptVec,:);
        eLengthVec=eLengthVec(isEKeptVec);
        nEdges=size(eMat,1);
        %% Build the first group of new edges (internal edges)
        e1NewMat=[...
            indVF12Vec,indVF23Vec;...     % 0
            indVF12Vec,indVF13Vec;...     % 1
            indVF23Vec,indVF13Vec];       % 2
        %% Build the first group of new faces (internal faces)
        f1NewMat=[indVF12Vec,indVF23Vec,indVF13Vec];
        %
        indNewEVec=transpose(nEdges+1:nEdges+nShrinkedFaces);
        %% Build face-to-edge map for the first group of new edges 
        %% and faces
        f2e1NewMat=repmat(indNewEVec,1,3)+...
            ones(nShrinkedFaces,1)*([0 2 1]*nShrinkedFaces);
        %% Build face-to-edge directions for the first group of new edges
        %% and faces
        f2e1IsDirNewMat=true(size(f2e1NewMat));
        %% Build the second group of new edges (edges on 
        %% the boundaries of the partitioned faces)
        indNewVert=transpose(nVerts+1:1:nVerts+nNewVerts);
        e2NewMat=[...
            indV1Vec,indNewVert;...
            indNewVert,indV2Vec];
        %% Build the second group of new faces (faces on the boundaries 
        %% of partitioned faces)
        f2NewMat=[...
            indVF12Vec,indVF13Vec,indVF1Vec;... 
            indVF23Vec,indVF12Vec,indVF2Vec;... 
            indVF13Vec,indVF23Vec,indVF3Vec];
        %% Build face-to-edge map for the second group of faces and edges
        indShiftEdgeVec=transpose(1:nShrinkedFaces);
        %
        dirMat=~xor(...
            kron([0 0;1 0;1 1]*nShrinkedFaces,ones(nShrinkedFaces,1)),...
            [f2eIsDirPartMat(:,[3 1]);...
            f2eIsDirPartMat(:,[1 2]);...
            f2eIsDirPartMat(:,[2 3])]);
        %
        %#1 - edges from the first group
        %#2 - shift indices for vertices correspond to shift idices for 
        %    edges
        f2e2NewMat=nEdges+...
            [repmat(indShiftEdgeVec,3,1)+kron([1;0;2]*nShrinkedFaces,...
            ones(nShrinkedFaces,1)),...%#1
            3*nShrinkedFaces-nVerts+dirMat*nNewVerts+... %#2
            [indVF13Vec indVF12Vec;indVF12Vec indVF23Vec;...
            indVF23Vec indVF13Vec]];
        %
        f2e2IsDirNewMat=logical(...
            [kron([1;0;0],ones(nShrinkedFaces,1)),dirMat]);
        %% Build the third group of new faces that correspond to edges 
        %% that have only 1 partitioned face of 2 adjacent faces
        % - this identifies edges in eMat
        indEBrokenButKeptVec=indEKeptVec(isE2PartVec&isEKeptVec);
        %but we still need to find indices of new vertices in the middles
        %of these edges
        indVMidEBrokenButKeptVec=indEPartVec(isE2PartVec&isEKeptVec);
        %faces have zero volume
        f3NewMat=[indVMidEBrokenButKeptVec,eMat(indEBrokenButKeptVec,:)];
        %
        ind12f2e3Vec=indVMidEBrokenButKeptVec-nVerts+3*nShrinkedFaces+nEdges;
        ind23f2e3Vec=ind12f2e3Vec+nNewVerts;
        %e2NewMat(indVMidEBrokenButKeptVec-nVerts,:)
        f2e3NewMat=[ind12f2e3Vec,ind23f2e3Vec,indEBrokenButKeptVec];
        f2e3IsDirNewMat=true(size(f2e3NewMat));
        %% Optionally udjust vertices
        if isAdjustFuncSpec
            vNewMat=fVertAdjustFunc(vNewMat);
        end
        %% Update edges, faces, vertices and f2eMap
        vMat=[vMat;vNewMat];
        eMat=[eMat;e1NewMat;e2NewMat];
        fMat=[fMat;f1NewMat;f2NewMat;f3NewMat];
        f2eMat=[f2eMat;f2e1NewMat;f2e2NewMat;f2e3NewMat];
        f2eIsDirMat=[f2eIsDirMat;f2e1IsDirNewMat;f2e2IsDirNewMat;...
            f2e3IsDirNewMat];
        %% Update edge length vec
        dMat=vMat(eMat(nEdges+1:1:end,1),:)-vMat(eMat(nEdges+1:1:end,2),:);
        eLengthVec=[eLengthVec;realsqrt(sum(dMat.*dMat,2))];
        %% Update number of entities
        nVerts=size(vMat,1);
        nEdges=size(eMat,1);
        nFaces=size(fMat,1);
        %% Collect stats
        if isStatCollected
            nVertVec=[nVertVec;nVerts];
            nEdgeVec=[nEdgeVec;nEdges];
            nFaceVec=[nFaceVec;nFaces];
        end
        iStep=iStep+1;        
    else
        break;
    end
end
if isStatCollected
    SStats.nSteps=iStep-1;
    SStats.nVertVec=nVertVec;
    SStats.nFaceVec=nFaceVec;
    SStats.nEdgeVec=nEdgeVec;
    SStats.nEdgesToShrinkVec=nEdgesToShrinkVec;
    SStats.maxEdgeLengthVec=maxEdgeLengthVec;
end