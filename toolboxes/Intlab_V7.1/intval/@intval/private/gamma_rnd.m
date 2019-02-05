function [y,ysup] = gamma_rnd(x,rnd)
% input x>=1 real column vector with gamma(xmax) overflow
% rnd  -1  y = lower bound for gamma(x)
%       1  y = upper bound for gamma(x)
%      []  [y,ysup] inclusion of gamma(x)
% rounding may be altered after leaving gamma_rnd
%

% written  06/19/13     S.M. Rump
%

  xmax = hex2num('406573fae561f647');       % gamma(x) overflow for x>xmax ~ 171.6
  
  y = x;
  if nargout==2
    ysup = x;
  end
  
  index = ( x>xmax );
  if any(index)
    y(index) = realmax;
    if nargout==2
      ysup(index) = inf;
    end
    index = ~index;
    if any(index)
      [y(index),ysup(index)] = gamma_rnd(x(index),rnd);
    end
    return
  end
  
  % 1 <= x <= xmax
  len = length(x);                          % number of elements
  col = floor(x) - 1;
  maxcol = max(col);
  factor = repmat(x,1,maxcol) - repmat(1:maxcol,len,1);
  factor(factor<1) = 1;
  if isempty(rnd)
    setround(-1)
    F = prod(factor,2);
    setround(1)
    Fsup = prod(factor,2);
  else
    setround(rnd)
    F = prod(factor,2);
  end
  
  x = x - col;                              % 1 <= x < 2, no rounding error
  gamma_data = getappdata(0,'INTLAB_INTVAL_GAMMADATA');
  AppGammaxs = gamma_data.AppGammaxs;       % approximate values gamma(xs)
  N = gamma_data.Nd(1);                     % number of grid points 1+i*h, h=1/N
  AppPoly = gamma_data.AppPoly;             % Approximating polynomials
  d = gamma_data.Nd(2);                     % degree of polynomials
  minerr = 0.560e-16;       % AppGammaxs-minerr <= gamma(xs) <= AppGammaxs+maxerr
  maxerr = 0.563e-16;
  ii = round(N*(x-1));
  xs = ii/N + 1;
  ii = ii + 1;
  delta = x - xs;
  if isempty(rnd) | ( rnd==-1 )
    setround(-1)
    y = AppPoly(ii,1);
    for k=2:d
      y = y.*delta + AppPoly(ii,k);
    end
    y = ( AppGammaxs(ii) + ( y.*delta - minerr ) ) .* F;
  end
  if isempty(rnd) | ( rnd==1 )
    setround(1)
    ysup = AppPoly(ii,1);
    for k=2:d
      ysup = ysup.*delta + AppPoly(ii,k);
    end
    ysup = ( AppGammaxs(ii) + ( ysup.*delta + maxerr ) );
    if rnd==1
      ysup = ysup .* F;
    else
      ysup = ysup .* Fsup;
    end
  end
  if rnd==1
    y = ysup;
  end
  