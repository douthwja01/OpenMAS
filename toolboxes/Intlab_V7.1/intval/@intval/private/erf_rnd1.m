function [y,ysup] = erf_rnd1(x,rnd,N)
% input x real non-negative column vector
% rnd  -1  y = lower bound for erf(x)
%       1  y = upper bound for erf(x)
%      []  [y,ysup] inclusion of erf(x)
% rounding may be altered after leaving erf_rnd1

% written  05/30/13     S.M. Rump
%

  % factorLB <= 2/sqrt(pi) <= factorUB
  INTLAB_STDFCTS_ERF = getappdata(0,'INTLAB_STDFCTS_ERF');
  factorLB = INTLAB_STDFCTS_ERF.TWO_SQRTPIINF;
  factorUB = INTLAB_STDFCTS_ERF.TWO_SQRTPISUP;

  phi = (1:(2*N+2));
  phi = (2*phi+1).*cumprod(phi);          % N<=9 to avoid rounding errors

  if isempty(rnd)
    Rnd = [-1 1];
  else
    Rnd = rnd;
  end
  for r=Rnd
    setround(r)
    x2n = (-x).*x;              % lower/upper bound for -x^2
    x4 = x.*x;                  % lower/upper bound for x^2
    x4 = x4.*x4;                % lower/upper bound for x^4
    if r==-1
      s = 0;
    else
      s = x4/phi(2*N+2);
    end
    for n=N:-1:1
      s = ( s + 1/phi(2*n) + x2n/phi(2*n+1) ) .* x4;
    end
    s = x2n/3 + s;
    if r==-1
      t = factorLB*x;
      y = t + t.*s;
    else                        % rnd==1 or []
      t = factorUB*x;
      if isempty(rnd)
        ysup = t + t.*s;
      else
        y = t + t.*s;
      end
    end
  end
