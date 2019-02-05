function r = isfinite(p)
%ISFINITE     Logical result:  polynomial p contains only finite components, y/n
%
%   r = isfinite(p)
%
%Result 1 iff all coefficients of p are finite
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = all(isfinite(p.c));
