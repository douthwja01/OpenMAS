function y = cos(x)
%COS          Implements  cos(x)  for intervals
%
%   y = cos(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, NaN input, sparse input,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 11/03/12     S.M. Rump  cos(0) (inspired by Sándor Kolumbán)
%

  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      y = intval(ones(size(x)));
      index = ~index;
      %VVVV  y(index) = cos(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,cos(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = cos(full(x));
    end
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
    y = ( exp(j*x) + exp(-j*x) ) / 2;
    if rndold
      setround(rndold)
    end
    return
  end

  if issparse(x.inf)                   % real case
    x.inf = full(x.inf);
    x.sup = full(x.sup);
  end

  y = x;

  % transform x.inf and x.sup mod pi/2
  [ xinfinf , xinfsup , Sinf ] = modpi2(x.inf(:));
  [ xsupinf , xsupsup , Ssup ] = modpi2(x.sup(:));
  Sinf = Sinf + 2;
  Sinf(Sinf>7) = Sinf(Sinf>7) - 8;
  Ssup = Ssup + 2;
  Ssup(Ssup>7) = Ssup(Ssup>7) - 8;

  [ yinf , ysup ] = sin_(x.inf(:),xinfinf,xinfsup,Sinf,  ...
    x.sup(:),xsupinf,xsupsup,Ssup );
  index = ~isfinite(x.inf(:)) | ~isfinite(x.sup(:));
  if any(index)
    yinf(index) = -1;
    ysup(index) = 1;
  end

  y.inf = reshape(yinf,size(x.inf));
  y.sup = reshape(ysup,size(x.inf));

  setround(0)                        % set rounding to nearest

  % take care of point interval zero
  index = ( ( x.inf==0 ) & ( x.sup==0 ) );
  if any(index(:))
    y.inf(index) = 1;
    y.sup(index) = 1;
  end
  
  index = isnan(x.inf);
  if any(index(:))
    y.inf(index) = NaN;
    y.sup(index) = NaN;
  end
  
  setround(rndold)
