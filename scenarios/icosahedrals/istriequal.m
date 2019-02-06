function [isPos,reportStr]=istriequal(v1Mat,f1Mat,v2Mat,f2Mat,maxTol)
% ISTRIEQUAL checks if two triangulations are equal. Triangulations
% resulting from permunations of edge directions and vertices in 
% faces are considered equal
% 
% Input:
%   regular:
%       v1Mat: double[n1Verts,3] - vertices of the first triangualtion
%       f1Mat: double[n1Faces,3] - faces of the first triangulation
%       v2Mat: double[n1Verts,3] - vertices of the second triangulation
%       f2Mat: double[n1Faces,3] - faces of the second triangulation
%       maxTol: double[1,1] - tolerances used for comparing vertex
%          coordinates
%
% Output:
%   isPos: logical[1,1] - specifies if result of comparison is true
%   reportStr: char[1,] - describes a reason of negative result
%
% $Author: Peter Gagarinov, PhD  <pgagarinov@gmail.com> $
% $Copyright: Peter Gagarinov, PhD, 
%            Moscow State University,
%            Faculty of Computational Mathematics and Computer Science,
%            System Analysis Department 2011-2016 $
%
n1Verts=size(v1Mat,1);
n2Verts=size(v2Mat,1);
isPos=n1Verts==n2Verts;
%
if isPos
    n1Faces=size(f1Mat,1);
    n2Faces=size(f2Mat,1);
    isPos=n1Faces==n2Faces;
    if isPos
        [v1Mat,~,indF1Vec]=sortrowstol(v1Mat,maxTol);
        [v2Mat,~,indF2Vec]=sortrowstol(v2Mat,maxTol);
        %
        f1Mat=indF1Vec(f1Mat);
        f2Mat=indF2Vec(f2Mat);
        realTol=max(max(abs(v1Mat-v2Mat)));
        isPos=realTol<=maxTol;
        %
        if isPos
            nF1Unique=size(unique(f1Mat,'rows'),1);
            nF2Unique=size(unique(f2Mat,'rows'),1);
            isPos=nF1Unique==nF2Unique;
            if isPos
                isPos=all(isface(v2Mat,f2Mat,f1Mat));
                if isPos
                    reportStr='';
                else
                    reportStr='faces are different';
                end
            else
                reportStr='numbers of unique faces are different';
            end
        else
            reportStr=sprintf(...
                'vertices are different, real tol=%f, exp tol=%f',...
                realTol,maxTol);
        end
    else
        reportStr='number of faces is different';
    end
else
    reportStr='numbers of vertices are different';
end