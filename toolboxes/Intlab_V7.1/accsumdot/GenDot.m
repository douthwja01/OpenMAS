function [x,y,d,cond] = GenDot(n,cond)
%GENDOT       Generation of extremely ill-conditioned dot products
%
%    [x,y,d,cond] = GenDot(n,cond)
%
%Generates (column) n-vectors x,y such that x'*y is approximately of condition cond.
%n should be at least 10.
%
%Input   n     length of vectors
%        cond  anticipated condition number
%
%Output  x,y   vectors
%        d     the value of x'*y rounded to double precision
%        cond  the true condition number of x'*y
%
%Implements Algorithm 6.1 in 
%   T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005
%

% written  11/08/03     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/01/08     S.M. Rump  comment
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/15/12     S.M. Rump  random generation
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  b = log2(cond);                       % b bits should cancel
  k = ceil(b/45)+1;
  n2 = round(n/2);
  x = zeros(n,1);
  y = x;
  
  e = round(rand(n2,1)*b/2);            % exponents between 0 and b/2
  e(1) = round(b/2)+1;                  % make sure exponents b/2 and
  e(end) = 0;                           %   0 actually occur
  x(1:n2) = randn(n2,1).*(2.^e);        % generate first half of vectors
  y(1:n2) = randn(n2,1).*(2.^e);
  
  e = round(linspace(b/2,0,n-n2));      % generate exponents for second half
  for i=n2+1:n
    x(i) = randn*2^e(i-n2);             % generate x(i), y(i) such that x(v)'*y(v)~2^e(i-n2)
    y(i) = (randn*2^e(i-n2)-DotK(x',y,k))/x(i);   % for v=1:i
  end
  
  index = randperm(n);                  % generate random permutation for x,y
  x = x(index);
  y = y(index);
  d = DotK(x',y,k);
  cond = 2*max(abs(x.*y))/abs(d);
    
  if rndold
    setround(rndold)
  end
