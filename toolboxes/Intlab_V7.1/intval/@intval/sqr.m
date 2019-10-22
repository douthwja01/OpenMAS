function y = sqr(x)
%SQR          Implements (elementwise)  sqr(x)  for intervals
%
%   y = sqr(x)
%

% written  12/06/98     S.M. Rump
% modified 12/24/99     S.M. Rump  sparse input, right bound zero
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,sqr(full(sx)),m,n);
    return
  end

  if x.complex
    y = x.*x;
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  y = x;

  setround(-1)
  y.inf = x.inf .* x.inf;

  index0 = ( x.inf<=0 ) & ( x.sup>=0 );
  if any(index0(:))
    y.inf(index0) = 0;
  end

  indexneg = ( x.sup<0 );
  if any(indexneg(:))
    y.inf(indexneg) = x.sup(indexneg) .* x.sup(indexneg);
  end

  setround(1)
  y.sup = x.sup .* x.sup;

  if any(index0(:))
    y.sup(index0) = max( y.sup(index0) , x.inf(index0) .* x.inf(index0) );
  end

  if any(indexneg(:))
    y.sup(indexneg) = x.inf(indexneg) .* x.inf(indexneg);
  end

  setround(rndold)
  