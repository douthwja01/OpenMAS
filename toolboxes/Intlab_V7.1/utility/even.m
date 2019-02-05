function res = even(n);
%EVEN         Boolean function: n even?
%
%   res = 1    n even
%   res = 0    n odd
%
%   res = even(n);
%

% written  10/29/97     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  res = ~mod(n,2);

