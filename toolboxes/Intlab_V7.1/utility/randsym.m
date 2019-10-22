function A = randsym(n,cnd)
%RANDSYM      Random symmetric matrix
%
%   res = randsym(n,cnd)
%
% cnd, if specified, is the approximate condition number of the generated matrix
%

% written  11/05/12     S.M. Rump
% modified 12/06/12     S.M. Rump  lower case
%

  if nargin==1
    A = randn(n);
  else
    s = exp( [ 0 rand(1,n-2) 1 ] * log(cnd) ) / cnd;
    A = randorth(n);
    A = A' * diag(s.*sign(randn(1,n))) * A;
  end

  % make sure A is symmetric
  A = tril(A)+tril(A,-1)';
