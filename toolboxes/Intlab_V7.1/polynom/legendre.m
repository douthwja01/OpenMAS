function p = legendre(n,var)
%LEGENDRE     Legendre polynomial
%
%   p = legendre(n,var)
%
%Computed by recursion
%  p(0,x) = 1
%  p(1,x) = x
%  p(n,x) = (2*n-1)/n * x * p(n-1,x) - (n-1)/n * p(n-2,x)    for n>1
%
%Parameter var for dependent variable optional (default 'x')
%

% written  07/30/02     S.M. Rump
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

  if nargin==1
    var = 'x';
  else
    if ~ischar(var)
      error('variable must be string')
    end
  end
  
  if n==0
    p = polynom(1,var);
  elseif n==1
    p = polynom([1 0],var);
  else
    t1 = [1 zeros(1,n)];
    t2 = [0 1 zeros(1,n-1)];
    for i=3:n+1
      t3 = (2*i-3)/(i-1)*[0 t2(1:n)] - (i-2)/(i-1)*t1;
      t1 = t2; 
      t2 = t3; 
    end
    p = polynom(fliplr(t3),var);
  end
  
  if rndold
    setround(rndold)
  end
