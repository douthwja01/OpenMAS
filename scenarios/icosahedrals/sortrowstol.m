function [resMat,indSortVec,indRevSortVec]=sortrowstol(inpMat,tol)
% SORTROWSTOL sorts rows of input numeric matrix in ascending order with a
% specified precision i.e. sorting [1 2;1+1e-14 1] with tol>=1-14 would put
% the first row on the second position while the built-in sortrows function
% would keep the order of rows unchanged, the functions only looks at the
% neighboring values and doesn't calculate pairwise distances (as in pdist
% function) for speed
%
% Input:
%   regular:
%       inpMat: numeric[nRows,nCols] - input matrix
%       tol: numeric[1,1] - tolerance used for sorting, values of inpMat
%          are considered to be equal if their difference is less or equal
%          than tol by an absolute value
%
% Output:
%   resMat: numeric[nRows,nCols] - output of matrix
%   indSortVec: double[nRows,1] - vector of indices such that
%       inpMat(indSortVec,:)==resMat
%   indRevSortVec: double[nRows,1] - vector of indices such that
%       resMat(indRevSortVec,:)==inpMat
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD,
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
% TODO: add support for clustering method based on
%   clusterdata([1;1+1e-14;2;2+1e-14;2-1e-14],'criterion',...
%       'distance','cutoff',1e-14)
%
if tol<0
    error('sortrowstol:wrongInput',...
        'tol is expected to be a positive number');
end
if ~(ismatrix(inpMat)&&isnumeric(inpMat))
    error('sortrowstol:wrongInput',...
        'input is expected to be a numeric matrix');
end
nCols=size(inpMat,2);
nRows=size(inpMat,1);
if nRows>0
    resMat=inpMat;%make a copy since inpMat is going to be changed
    %
    %cluster values in each column
    for iCol=1:nCols
        [colVec,indColSortVec]=sort(inpMat(:,iCol));
        colDiffVec=diff(colVec);
        isLessVec=abs(colDiffVec)<=tol;
        colDiffVec(isLessVec)=0;
        colVec=cumsum([colVec(1);colDiffVec],1);
        [~,indColRevSortVec]=sort(indColSortVec);
        inpMat(:,iCol)=colVec(indColRevSortVec);
    end
    %sort the original matrix copied to resMat
    [~,indSortVec]=sortrows(inpMat);
    resMat=resMat(indSortVec,:);
else
    resMat=inpMat;
    indSortVec=double.empty(0,1);
end
if nargout>2
    [~,indRevSortVec]=sort(indSortVec);
end