function c = abs(a)
%ABS          Interval absolute value
%
%  c = abs(a);
%
%Result c is the real interval of { abs(alpha) | alpha in a }
%

% written  12/06/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 07/08/02     S.M. Rump  changed from iabs-> abs
% modified 04/04/04     S.M. Rump  redesign, set round to nearest for safety
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 06/26/06     S.M. Rump  zero bound corrected (thanks to Matthew Benedict)
% modified 09/26/06     S.M. Rump  zero bound corrected
% modified 10/03/09     S.M. Rump  sparse input
% modified 08/17/12     S.M. Rump  treatment of NaN (thanks to Nikolay Petrov 
%                                    for pointing to this)
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  c = a;
  c.complex = 0;
  
  if a.complex
    
    setround(-1)
    if isequal(a.rad,0)            % take care of huge arrays
      c.inf = abs(a.mid);
      setround(1)
      c.sup = abs(a.mid);
    else
      c.inf = abs(a.mid) - a.rad;
      c.inf(c.inf<0) = 0;
      setround(1)
      c.sup = abs(a.mid) + a.rad;
    end
    
    if issparse(a.mid)    % take care of huge arrays
      [i,j] = find(isnan(a.mid+a.rad));
      [m n] = size(a.mid);
      c.inf = c.inf + sparse(i,j,NaN,m,n);
      c.sup = c.sup + sparse(i,j,NaN,m,n);
    else
      index = find(isnan(a.mid));
      if any(index(:))
        c.inf(index) = NaN;
        c.sup(index) = NaN;
      end
    end

  else
    
    c.inf = min(abs(a.inf),abs(a.sup));
    c.sup = max(abs(a.inf),abs(a.sup));
    
    index = find( ( a.inf<0 ) & ( a.sup>0 ) );
    if any(index(:))
      if issparse(c.inf)
        c.inf(index) = complex(0,1);    % here come the dirty tricks
        c.inf = real(c.inf);
      else
        c.inf(index) = 0;   % careful: terribly slow for large sparse input
      end
    end
    
    if issparse(a.inf)    % take care of huge arrays
      [i,j] = find(isnan(a.inf));
      [m n] = size(a.inf);
      c.inf = c.inf + sparse(i,j,NaN,m,n);
      c.sup = c.sup + sparse(i,j,NaN,m,n);
    else
      index = isnan(c.inf);
      if any(index(:))
        c.inf(index) = NaN;
        c.sup(index) = NaN;
      end
    end
    
  end
  
  c.mid = [];
  c.rad = [];
  
  setround(rndold)
