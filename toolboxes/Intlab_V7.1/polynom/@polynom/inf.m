function p = inf(p);
%INF          Polynomial of lower bounds for (interval) polynomial (same as p.inf)
%
%   r = inf(p)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = inf_(p.c);
  p = normalize(p);
