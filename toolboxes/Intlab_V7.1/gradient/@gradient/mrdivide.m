function r = mrdivide(a,b)
%MRDIVIDE     Implements  a / b  for gradient
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    complete redesign
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if prod(size(b))~=1
    error('Gradient division only for scalar denominator')
  end

  N = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  if ~isa(a,'gradient')          % non-gradient / gradient scalar
    r.x = a / b.x;
    n = prod(size(a));
    if n==1                     % non-gradient scalar / gradient scalar
      D = -a / sqr(b.x);
      r.dx = b.dx * D;
    else                        % non-gradient array / gradient scalar
      D = 1 / sqr(b.x);
      if issparse(b.dx)
        ax = sparse(-a(:));
      else
        ax = -a(:);
      end
      r.dx = ax * ( b.dx * D );
    end
  elseif ~isa(b,'gradient')      % gradient / non-gradient scalar    
    r.x = a.x / b;
    r.dx = a.dx / b;
  else                          % gradient / gradient scalar
    r.x = a.x / b.x;
    D = 1 / sqr(b.x);
    n = prod(size(a.x));
    if n==1                     % gradient scalar / gradient scalar
      Num = a.dx * b.x - a.x * b.dx;
      r.dx = Num * D;
    else                        % gradient array / gradient scalar
      if issparse(a.dx)
        ax = sparse(a.x(:));
      else
        ax = a.x(:);
      end
      Num = a.dx * b.x - ax * b.dx;
      r.dx = Num * D;
    end
  end

  r = class(r,'gradient');
  
  if rndold
    setround(rndold)
  end
