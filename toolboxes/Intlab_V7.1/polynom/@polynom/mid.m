function p = mid(p);
%MID          Midpoint of (interval) polynomial (same as p.mid)
%
%   r = mid(p)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = mid(p.c);
  p = normalize(p);
