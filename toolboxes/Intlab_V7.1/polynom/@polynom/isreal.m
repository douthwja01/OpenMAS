function r = isreal(p)
%ISREAL       Boolean function
%
%   r = isreal(p);
%
%Result 1 iff all coefficients of p are real
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = isreal(p.c);
