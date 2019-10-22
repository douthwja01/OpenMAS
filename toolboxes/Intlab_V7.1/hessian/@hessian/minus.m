function r = minus(a,b)
%MINUS        Hessian addition  a - b
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/23/06     S.M. Rump  scalar+array (thanks to Sébastien Loisel)
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(a,'hessian')          % non-hessian plus hessian
    r.x = a - b.x;
    if prod(size(b.x))>1        % b is not scalar
      r.dx = - b.dx;
      r.hx = - b.hx;
    else                        % b is scalar, a may be array
      M = prod(size(r.x));
      r.dx = - b.dx(:,ones(1,M));
      r.hx = - b.hx(:,ones(1,M));
    end
    if isa(a,'intval') & isa(b.x,'double')
      r.dx = intval(r.dx);
      r.hx = intval(r.hx);
    end
  elseif ~isa(b,'hessian')      % hessian plus non-hessian
    r.x = a.x - b;
    if prod(size(a.x))>1        % a is not scalar
      r.dx = a.dx;
      r.hx = a.hx;
    else                        % a is scalar, b may be array
      M = prod(size(r.x));
      r.dx = a.dx(:,ones(1,M));
      r.hx = a.hx(:,ones(1,M));
    end
    if isa(b,'intval') & isa(a.x,'double')
      r.dx = intval(r.dx);
      r.hx = intval(r.hx);
    end
  else                          % both input a and b are Hessian
    M = prod(size(a));
    N = prod(size(b));
    r.x = a.x - b.x;
    if ( M==1 ) & ( N~=1 )      % scalar Hessian - array Hessian
      r.dx = a.dx(:,ones(1,N)) - b.dx;
      r.hx = a.hx(:,ones(1,N)) - b.hx;
    elseif ( M~=1 ) & ( N==1 )  % array Hessian - scalar Hessian
      r.dx = a.dx - b.dx(:,ones(1,M));
      r.hx = a.hx - b.hx(:,ones(1,M));
    else                        % addition of both scalar or both array Hessians
      r.dx = a.dx - b.dx;
      r.hx = a.hx - b.hx;
    end
  end

  r = class(r,'hessian');
  
  if rndold
    setround(rndold)
  end
