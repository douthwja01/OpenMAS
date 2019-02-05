function r = isinf(p)
%ISINF        Logical result:  polynomial p contains inf components, y/n
%
%   r = isinf(p)
%
%Result 1 iff at least one coefficient of p contains some inf
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = any(isinf(p.c));
