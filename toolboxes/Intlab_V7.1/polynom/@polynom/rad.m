function p = rad(p);
%RAD          Radius of (interval) polynomial (same as p.rad)
%
%   r = rad(p)
%

% written  10/04/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = rad(p.c);
  p = normalize(p);
