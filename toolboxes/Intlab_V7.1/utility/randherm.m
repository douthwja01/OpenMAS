function A = randherm(n,cnd)
%RANDHERM     Random Hermitian matrix
%
%   res = randherm(n,cnd)
%
% cnd, if specified, is the approximate condition number of the generated matrix
%

% written  11/05/12     S.M. Rump
% modified 12/06/12     S.M. Rump  lower case, typo
%

  if nargin==1
    A = randn(n) + sqrt(-1)*randn(n);
  else
    s = exp( [ 0 rand(1,n-2) 1 ] * log(cnd) ) / cnd;
    A = orth(randn(n)+sqrt(-1)*randn(n));
    A = A' * diag(s.*sign(randn(1,n))) * A;
  end

  % make sure A is Hermitian
  A = A - sqrt(-1)*diag(imag(diag(A)));   % make sure diagonal is real
  A = tril(A)+tril(A,-1)';
