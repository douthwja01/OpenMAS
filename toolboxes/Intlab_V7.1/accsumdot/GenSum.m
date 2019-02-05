function [x,d,cond] = GenSum(n,cond)
%GENSUM       Generation of extremely ill-conditioned sums
%
%    [x,d,cond] = GenSum(n,cond)
%
%Generates (column) n-vector x such that sum(x) is approximately of condition cond.
%n should be at least 10.
%
%Input   n     length of vector
%        cond  anticipated condition number
%
%Output  x     vector
%        d     the value of sum(x) rounded to double precision
%        cond  the true condition number of sum(x)
%
%Adapted from Algorithm 6.1 in 
%   T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005
%

% written  12/10/05     S.M. Rump
% modified 09/01/08     S.M. Rump  comment
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/15/12     S.M. Rump  random generation
% modified 12/05/12     S.M. Rump  typo
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
  
  e = round(rand(n2,1)*b);              % exponents between 0 and b
  e(1) = round(b/2)+1;                  % make sure exponents b/2 and
  e(end) = 0;                           %   0 actually occur
  x(1:n2) = randn(n2,1).*(2.^e);        % generate first half of vector
  
  e = round(linspace(b,0,n-n2));        % generate exponents for second half
  for i=n2+1:n
    % generate x(i) such that sum(x(v))~2^e(i-n2) for v=1:i
    x(i) = (randn*2^e(i-n2)-SumK(x,k));
  end
  
  index = randperm(n);                  % generate random permutation for x,y
  x = x(index);
  d = SumK(x,k);
  cond = sum(abs(x))/abs(d);
    
  if rndold
    setround(rndold)
  end
