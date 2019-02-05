function p = sup(p);
%SUP          Polynomial of upper bounds for (interval) polynomial (same as p.sup)
%
%   r = sup(p)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = sup(p.c);
  p = normalize(p);
