function [y,ysup] = erfc_rnd3(x,rnd,N)
% input x real non-negative column vector
%   xmax = hex2num('403b369a6244e684');   % ~27.2: erfc(x)<subrealmin fuer x>=xmax
% rnd  -1  y = lower bound for erfc(x)
%       1  y = upper bound for erfc(x)
%      []  [y,ysup] inclusion of erfc(x)
% rounding may be altered after leaving erfc_rnd3
% 

% written  05/30/13     S.M. Rump
%

  % factorLB <= 2/sqrt(pi) <= factorUB
  INTLAB_STDFCTS_ERF = getappdata(0,'INTLAB_STDFCTS_ERF');
  factorLB = INTLAB_STDFCTS_ERF.ONE_SQRTPIINF;
  factorUB = INTLAB_STDFCTS_ERF.ONE_SQRTPISUP;  % ~ 0.56
  
  setround(0)
  factor = hex2num('41a0000002000000');     % 2^27+1
  b = factor*x;                             % split x = a+b
  a = b - ( b - x );
  b = x - a;

  K = ceil(N/2);
  if isempty(rnd)
    Rnd = [-1 1];
  else
    Rnd = rnd;
  end
  for r=Rnd
    setround(r)
    x_2 = - ( ((-0.5)./x)./x );     % upper/lower bound of 1/2x^2
    x_4 = 0.5./x./x;                % lower/upper bound of 1/2x^2
    x_4 = x_4.*x_4;                 % lower/upper bound of 1/4x^4
    if r==-1
      s = 0;
    else
      s = ( (4*K+3)*(4*K+1) ) * x_4;
    end
    for k=K:-1:1
      s = ( 1 + ( s + (-4*k-1)*x_2 ) ) .* ( ((4*k-1)*(4*k-3))*x_4 );
    end
    % a^2, a*b and b^2 computed w/o rounding error
    phi = exp_rnd(-a.*a,r).*exp_rnd((-2*a.*b)-(b.*b),r);  % bound for exp(-x^2)
    if r==-1
      y = phi .* ( factorLB * ( ( 1 + ( s - x_2 ) ) ./ x ) );
    else
      if isempty(rnd)
        ysup = phi .* ( factorUB * ( ( 1 + ( s - x_2 ) ) ./ x ) );
      else
        y = phi .* ( factorUB * ( ( 1 + ( s - x_2 ) ) ./ x ) );
      end
    end
  end
