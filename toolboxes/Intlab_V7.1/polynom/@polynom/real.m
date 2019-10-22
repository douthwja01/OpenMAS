function p = real(p);
%REAL         Real part of (interval) polynomial
%
%   r = real(p)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = real(p.c);
  p = normalize(p);
