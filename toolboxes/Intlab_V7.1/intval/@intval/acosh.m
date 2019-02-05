function y = acosh(x)
%ACOSH        Implements  acosh(x)  for intervals
%
%   y = acosh(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, array input, sparse input
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 01/20/03     S.M. Rump  Matlab sqrt fixed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and
%                                     some improvements, extreme values for 
%                                     approximate part
% modified 09/06/07     S.M. Rump  approximate std fcts removed, exceptional arguments
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/20/08     S.M. Rump  check for zero omitted
% modified 10/18/08     S.M. Rump  StdFctsException ignore/NaN
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
    if issparse(x.mid)
      x.mid = full(x.mid);
      x.rad = full(x.rad);
    end
    if rndold
      setround(rndold)
    end
    y = 2 * log( sqrt((x+1)/2) + sqrt((x-1)/2) );  
    return
  end

  if issparse(x.inf)
    x.inf = full(x.inf);
    x.sup = full(x.sup);
  end
  % input x real and full
  % real range of definition:  [1,inf]
  INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
  index = ( x.inf<1 );                  % (partially) exceptional indices
  if any(index(:))                      % handle input out-of-range
    if INTLAB_STDFCTS_EXCPTN<=1  % out-of-range input handled as complex
      if INTLAB_STDFCTS_EXCPTN==1
        warning('ACOSH: Real interval input out of range changed to be complex')
      end
      y = x;
      %VVVV  y(index) = acosh(cintval(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acosh(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~index;
      if any(index(:))
        %VVVV  y(index) = acosh(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acosh(subsref(x,s)));
        %AAAA  Matlab bug fix
      end
      if rndold
        setround(rndold)
      end
      return
    end
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',1);
    if INTLAB_STDFCTS_EXCPTN==3    % ignore input out of range (ignore-mode)
      x.inf(index) = 1;                   % completely exceptional indices treated below
      index1 = index & ( x.sup<1);        % completely exceptional indices
    end
  else
    index1 = [];                          % make sure indexneg is not undefined
  end

  % input x real and full
  y = x;
  wng = warning;
  warning off
  
  % treat non-exceptional arguments
  y.inf = acoshrnd(x.inf,-1);
  y.sup = acoshrnd(x.sup,1);

  if INTLAB_STDFCTS_EXCPTN==3      % ignore input out of range (ignore-mode)
    if ~isempty(find(index1))             % completely exceptional arguments to NaN
      y.inf(index1) = NaN;
      y.sup(index1) = NaN;
    end
  else                                    % any input out of range to NaN (NaN-mode)
    if ~isempty(find(index))              % exceptional arguments to NaN
      y.inf(index) = NaN;
      y.sup(index) = NaN;
    end
  end

  setround(rndold)
  warning(wng)

  

function y = acoshrnd(x,rnd)
% rigorous acosh for nonnegative x with rounding corresponding to rnd

  y = x;

  % huge arguments: use proper scaling
  index = ( x>1e10 );
  if any(index(:))
    setround(rnd)
    y(index) = log_rnd( 2*x(index) + eps*(rnd-1) , rnd );
  end
  Index = ~index;

  % large arguments: use acosh(x) = log( x + sqrt(x^2-1) )
  index = Index & ( x>1.25 );
  if any(index(:))
    X = x(index);
    setround(rnd)
    y(index) = log_rnd( X + sqrt_rnd(sqr(X)-1,rnd) , rnd );
  end
  Index = Index & ( ~index );

  % small arguments: use acosh(x) = log( x + sqrt(E*(2+E)) )  for  x = 1+E
  index = Index & ( ~index );          % 1 <= x <= 1.25
  if any(index(:))
    X = x(index);
    E = X-1;              % exactly representable because X close to 1
    setround(rnd)
    E = E + sqrt_rnd( E.*(2+E) , rnd );
    y(index) = log_1(E,rnd);     % 0 <= E <= 1
  end
