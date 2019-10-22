function c = sum(a,dim)
%SUM          Implements  sum(a,dim)  for slopes
%
%   u = sum(a,dim)
%
%parameter dim optional, functionality as Matlab function sum
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==1
    if a.size(1)==1
      dim = 2;
    else
      dim = 1;
    end
  end

  if dim==1
    c = ones(1,a.size(1))*a;
  else
    c = a*ones(a.size(2),1);
  end
  
  setround(rndold)
