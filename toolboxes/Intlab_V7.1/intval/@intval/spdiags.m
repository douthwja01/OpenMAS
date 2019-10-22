function [C,d] = spdiags(A,d,B,n)
%SPDIAGS      Implements  spdiags  for intervals
%
% functionality as Matlab function spdiags for matrices
%

% written  07/08/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 02/16/11     S.M. Rump  call with 3 arguments (thanks to Dimitar Dimitrov)
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==1              % C = spdiags(A)  or  [C,d] = spdiags(A), A must be intval
    C.complex = A.complex;
    if C.complex            % A complex intval
      C.inf = [];
      C.sup = [];
      [C.mid,d] = spdiags(A.mid);
      if isequal(A.rad,0)
        C.rad = 0;
      else
        [C.rad,e] = spdiags(A.rad);
        if ~isequal(d,e)
          d = union(d,e);
          C.mid = spdiags(A.mid,d);
          C.rad = spdiags(A.rad,d);
        end
      end
    else                    % A real intval
      [C.inf,d] = spdiags(A.inf);
      [C.sup,e] = spdiags(A.sup);
      if ~isequal(d,e)
        d = union(d,e);
        C.inf = spdiags(A.inf,d);
        C.sup = spdiags(A.sup,d);
      end
      C.mid = [];
      C.rad = [];
    end
  elseif nargin==2        % C = spdiags(A,d)
    if isa(d,'intval')
      error('invalid call of intval/spdiags')
    end
    C.complex = A.complex;
    if C.complex
      C.inf = [];
      C.sup = [];
      C.mid = spdiags(A.mid,d);
      if isequal(A.rad,0)
        C.rad = 0;
      else
        C.rad = spdiags(A.rad,d);
      end
    else
      C.inf = spdiags(A.inf,d);
      C.sup = spdiags(A.sup,d);
      C.mid = [];
      C.rad = [];
    end
  elseif nargin==3        % C = spdiags(A,d,B)
    if isa(d,'intval')
      error('invalid call of intval/spdiags')
    end
    A = intval(A);
    B = intval(B);
    if A.complex | B.complex
      C.complex = 1;
      A = cintval(A);
      B = cintval(B);
    else
      C.complex = 0;
    end
    if C.complex
      C.inf = [];
      C.sup = [];
      C.mid = spdiags(A.mid,d,B.mid);
      if isequal(A.rad,0)
        C.rad = 0;
      else
        C.rad = spdiags(A.rad,d,B.rad);
      end
    else
      C.inf = spdiags(A.inf,d,B.inf);
      C.sup = spdiags(A.sup,d,B.sup);
      C.mid = [];
      C.rad = [];
    end
  elseif nargin==4        % C = spdiags(A,d,m,n), parameter m is B
    if isa(d,'intval')
      error('invalid call of intval/spdiags')
    end
    C.complex = A.complex;
    if C.complex
      C.inf = [];
      C.sup = [];
      C.mid = spdiags(A.mid,d,B,n);
      if isequal(A.rad,0)
        C.rad = 0;
      else
        C.rad = spdiags(A.rad,d,B,n);
      end
    else
      C.inf = spdiags(A.inf,d,B,n);
      C.sup = spdiags(A.sup,d,B,n);
      C.mid = [];
      C.rad = [];
    end
  end
  
  C = class(C,'intval');
    
  if rndold
    setround(rndold)
  end
