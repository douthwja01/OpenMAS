function r = inv(a)
%INV          Matrix inverse for interval square matrices
%
%  r = inv(a)
%
%Calls verifylss(a,speye(size(a)))
%

% written  10/16/98     S.M. Rump
%

  n = dim(a);                  % checks a to be square
  r = verifylss(a,speye(n));
