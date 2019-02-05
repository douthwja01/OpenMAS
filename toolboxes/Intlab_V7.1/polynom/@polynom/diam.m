function p = diam(p);
%DIAM         Diameter of (interval) polynomial
%
%   r = diam(p)
%

% written  10/04/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = diam(p.c);
  p = normalize(p);
