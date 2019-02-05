function c = sum(a,dim)
%SUM          Implements  sum(a,dim)  for intervals
%
%   c = sum(a,dim)
%
% parameter dim optional
% functionality as Matlab function sum for matrices
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
  end

  c.complex = a.complex;

  if a.complex
    if nargin==1,
      if size(a.mid,1)==1
        dim = 2;
      else
        dim = 1;
      end
    end
    c.inf = [];
    c.sup = [];
    setround(-1)
    c1 = sum(a.mid,dim);
    setround(1)
    c2 = sum(a.mid,dim);
    c.mid = c1 + 0.5*(c2-c1);
    if isequal(a.rad,0)
      c.rad = abs(c.mid-c1);
    else
      c.rad = abs(c.mid-c1) + sum(a.rad,dim);
    end
  else
    if nargin==1,
      if size(a.inf,1)==1
        dim = 2;
      else
        dim = 1;
      end
    end
    setround(-1)
    c.inf = sum(a.inf,dim);
    setround(1)
    c.sup = sum(a.sup,dim);
    c.mid = [];
    c.rad = [];
  end

  c = class(c,'intval');

  setround(rndold)                      % set rounding to previous value
