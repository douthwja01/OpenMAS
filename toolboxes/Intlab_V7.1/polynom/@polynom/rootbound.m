function r = rootbound(p)
%ROOTBOUND    Root bound for univariate (interval) polynomial
%
%   r = rootbound(p)
%

% written  09/14/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 05/31/09     S.M. Rump  interval polynomials
%

  if size(p.e,2)>1
    error('rootbound only for univariate polynomials')
  end

  n = degree(p);
  lc = mig(p.c(1));
  if lc==0              % no root bound possible, leading coefficient contains zero
    r = NaN;
    return
  end
  if n==0               % constant polynomial, must be nonzero
    r = 0;
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  pabs = mag(p.c);
  pabs(1) = lc;
  rb1 = 1 + max( pabs/lc );
  v = 1:n; 
  rb2 = 2 * max( ( pabs(v+1)/lc ) .^ ( 1./v ) );
  r = min(rb1,rb2);
  q = [lc -pabs(2:n+1)];
  if isa(p.c,'intval')
    rold = 0;
    Q = polynom(intval(q));
    Qs = pderiv(Q);
    while abs(r-rold) > 1e-4*abs(r)
      rold = r;
      corr = polyval(Q,r)/polyval(Qs,r);      % takes care of large r
      if ~isnan(corr)
        r = sup(r - corr);
      end
    end
  else
    qs = ( n:-1:1 ) .* q(1:n);                % derivative q'
    rold = 0;
    while abs(r-rold) > 1e-4*abs(r)
      rold = r;
      corr = polyval(q,r)/polyval(qs,r);      % takes care of large r
      if ~isnan(corr)
        r = r - corr;
      end
    end
  end

  if rndold
    setround(rndold)
  end
