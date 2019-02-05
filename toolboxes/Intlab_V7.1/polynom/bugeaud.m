function P = bugeaud(n,a,k)
%BUGEAUD      Bugeaud/Mignotte polynomial with very small root separation
%
%   P = bugeaud(n,a,k)
%
%defined by   P(X) = (X^n-aX+1)^k - 2X^(nk-k)*(aX-1)^k  for n>=3, a>=10 and k>=2.
%Default is a=10 and k=2.
%The polynomial has k roots in the circle with center 1/a+a^(-n-1) and radius 2a^(-2n).
%See Y. Bugeaud, M. Mignotte: On the distance between roots of integer polynomials,
%  Proceedings of the Edinburgh Mathematical Society 47: 553-556, 2004.
%

% written  04/13/08     S.M. Rump
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                   % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  if nargin<3                   % set defaults
    k = 2;
    if nargin<2
      a = 10;
    end
  end
  
  X = polynom([1 0]);
  P = (X^n - a*X + 1)^k - 2*X^(n*k-k)*(a*X-1)^k;
  
  if rndold
    setround(rndold)
  end
