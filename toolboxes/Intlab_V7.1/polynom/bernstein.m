function B = bernstein(n,i,var)
%BERNSTEIN    The i-th Bernstein polynomial B of degree n
%
%   B = bernstein(n,i,var)
%
%defined by   binom(n,i) * x^i * (1-x)^(n-i) .
%Input n and i may be vectors (of equal length with 0 <= i <= n).
%Parameter var for depending variable of Bernstein polynomial is optional, 
%  default is 'x' for univariate and 'x1',...,'xn' for multivariate polynomials.
%

% written  10/25/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  k = length(n);            % number of variables
  if nargin==2
    if k==1
      var = 'x';
    else
      for j=1:k
        var{j} = ['x' int2str(j)];
      end
    end
  end
  
  if k==1                   % univariate case
    
    j = (n-i):-1:0;
    f = -ones(1,n-i+1);
    f(end:-2:1) = 1;
    B.c = binom(n,i) * [f.*binom(n-i,j) zeros(1,i)];
    B = polynom(B.c,var);
    
  else                      % multivariate case
    
    B = bernstein(n(1),i(1),var{1});
    for j=2:k
      B = B * bernstein(n(j),i(j),var{j});
    end
    
  end
  
  if rndold
    setround(rndold)
  end
