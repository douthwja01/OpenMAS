function r = isnan(p)
%ISNAN        Logical result:  polynomial p contains NaN components, y/n
%
%   r = isnan(p)
%
%Result 1 iff at least one coefficient of p contains some NaN
%

% written  09/16/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = any(isnan(p.c));
