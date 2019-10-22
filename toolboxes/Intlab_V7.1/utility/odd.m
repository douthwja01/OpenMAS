function res = odd(n);
%ODD          Boolean function: n odd?
%
%   res = 1    n odd
%   res = 0    n even
%
%   res = odd(n);
%

% written  10/29/97     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  res = mod(n,2);

