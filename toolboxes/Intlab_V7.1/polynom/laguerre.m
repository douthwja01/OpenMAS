function p = laguerre(n,a,var)
%LAGUERRE     Laguerre polynomial
%
%   p = laguerre(n)
%   p = laguerre(n,var)
%   p = laguerre(n,a,var)
%
%Computed by recursion
%  p(0,a,x) = 1
%  p(1,a,x) = -x+1+a
%  p(n,a,x) = (2*n+a-1-x)/n*p(n-1,a,x) - (n+a-1)/n*p(n-2,a,x)  for n>1
%
%Input parameter "a" real number greater than -1, default a=0. 
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
    a = 0;
    var = 'x';
  elseif nargin==2
    if ischar(a)
      var = a;
      a = 0;
    else
      var = 'x';
    end
  else
    if ~ischar(var)
      error('variable must be string')
    end
  end
  
  if a<=-1
    error('input parameter a must be greater than -0.5')
  end
  
  if n==0
    p = polynom(1,var);
  elseif n==1
    p = polynom([-1 1+a],var);
  else
    t1 = [1 zeros(1,n)];
    t2 = [1+a -1 zeros(1,n-1)];
    for i=3:n+1
      t3 = (2*(i-1)+a-1)/(i-1)*t2 - [0 t2(1:n)]/(i-1) - (i+a-2)/(i-1)*t1; 
      t1 = t2; 
      t2 = t3; 
    end
    p = polynom(fliplr(t3),var);
  end
  
  if rndold
    setround(rndold)
  end
