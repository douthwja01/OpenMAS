function r = randint(n,k,l);
%RANDINT      Random integer (matrix) uniformly distributed
%
%n>0       uniformly in [1,n]
%otherwise uniformly in [-n,n]
%
%   r = randint(n);       random integer in range
%
%   randint(n,k)          k x k  matrix with entries randint(n)
%
%   randint(n,k,l)        k x l  matrix with entries randint(n)
%

% written  02/08/93     S.M. Rump
% modified 10/29/93     S.M. Rump
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

  if n>0
    if nargin==1
      r = floor(n*rand)+1;
    elseif nargin==2
      r = floor(n*rand(k))+1;
    elseif nargin==3
      r = floor(n*rand(k,l))+1;
    end
  else
    n=-n;
    if nargin==1
      r = floor((2*n+1)*rand)-n;
    elseif nargin==2
      r = floor((2*n+1)*rand(k))-n;
    elseif nargin==3
      r = floor((2*n+1)*rand(k,l))-n;
    end
  end
  
  if rndold
    setround(rndold)
  end
