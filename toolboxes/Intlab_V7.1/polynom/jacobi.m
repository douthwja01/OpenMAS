function p = jacobi(n,a,b,var)
%JACOBI       Jacobi polynomial
%
%   p = jacobi(n,a,b,var)
%
%Computed by recursion
%  p(0,a,b,x) = 1
%  p(1,a,b,x) = (a-b)/2 + (1+(a+b)/2)*x
%  p(n,a,b,x) = (2*n+a+b-1) * (a^2 - b^2 + (2*n+a+b-2)*(2*n+a+b)*x) / (2*n*(n+a+b)*(2*n+a+b-2)) * p(n-1,a,b,x)
%                 - 2*(n+a-1)*(n+b-1)*(2*n+a+b) / (2*n*(n+a+b)*(2*n+a+b-2)) * p(n-2,a,b,x)    for n>1
%
%Input parameters "a" and "b" real number greater than -1. 
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

  if nargin==3
    var = 'x';
  else
    if ~ischar(var)
      error('variable must be string')
    end
  end
  
  if ( a<=-1 ) | ( b<=-1 )
    error('input parameters a and b must be greater than -1')
  end
  
  if n==0
    p = polynom(1,var);
  elseif n==1
    p = polynom([1+(a+b)/2 (a-b)/2],var);
  else
    t1 = [1 zeros(1,n)];
    t2 = [(a-b)/2 1+(a+b)/2 zeros(1,n-1)];
    for i=3:n+1
      f1 = 2*i+a+b-3;
      f2 = 2*(i-1)*(i-1+a+b)*(2*i+a+b-4);
      t3 = f1*(a^2-b^2)/f2*t2 + f1*(2*i+a+b-4)*(2*i-2+a+b)/f2*[0 t2(1:n)] - ...
             2*(i+a-2)*(i+b-2)*(2*i-2+a+b)/f2*t1; 
      t1 = t2; 
      t2 = t3; 
    end
    p = polynom(fliplr(t3),var);
  end
  
  if rndold
    setround(rndold)
  end
