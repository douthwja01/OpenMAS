function isFaceVec=isface(vMat,fMat,fToCheckMat)
% ISFACE checks if the specified faces from fToCheckMat belong to the 
% given triangulation fMat of vertices from vMat
%
% Input:
%   regular:
%       vMat: double[nVerts,3] - vertex coordinates
%       fMat: double[nFaces,3] - face definitions based on vertex numbers
%       fToCheckMat: double[nCheckFaces,3] - definitions of faces which
%          beloning is checked
%
% Output:
%   isFaceVec: logical[nCheckFaces,1] - contains true for the corresponding
%       face from fToCheckMat if the face belongs to the triangulation and
%       false otherwise
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
tr=triangulation(fMat,vMat);
f12NumCVec=tr.edgeAttachments(fToCheckMat(:,[1,2]));
f23NumCVec=tr.edgeAttachments(fToCheckMat(:,[2,3]));
f13NumCVec=tr.edgeAttachments(fToCheckMat(:,[1,3]));
%
f12NumMat=fnumcvec2mat(f12NumCVec,[-121 -122]);
f23NumMat=fnumcvec2mat(f23NumCVec,[-231 -232]);
f13NumMat=fnumcvec2mat(f13NumCVec,[-131 -132]);
%
isFaceVec=intersectmat(f12NumMat,f23NumMat,f13NumMat);
end
%
function isPosRes=intersectmat(oneMat,twoMat,threeMat)
combMat=combvec([0,1],[0,1],[0,1]).';
nCombs=size(combMat,1);
isPosRes=false(size(oneMat,1),1);
for iComb=1:nCombs
    isPosRes=isPosRes|cmpvec(combMat(iComb,:));
end
    function isPos=cmpvec(isFlipVec)
        isPos=any(...
            circshift(oneMat,[0 isFlipVec(1)])==circshift(twoMat,...
            [0 isFlipVec(2)])&...
            circshift(twoMat,[0 isFlipVec(2)])==circshift(threeMat,...
            [0 isFlipVec(3)]),2);
    end
end
%
function fNumMat=fnumcvec2mat(fNumCVec,defVec)
nElems=size(fNumCVec,1);
lVec=cellfun('length',fNumCVec);
fNumMat=repmat(defVec,nElems,1);
isOneElemVec=lVec==1;
fNumMat(isOneElemVec,1)=[fNumCVec{isOneElemVec}].';
isTwoElemVec=lVec==2;
fNumMat(isTwoElemVec,:)=reshape([fNumCVec{isTwoElemVec}],2,[]).';
end