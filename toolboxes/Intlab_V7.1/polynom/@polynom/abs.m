function p = abs(p);
%ABS          Polynomial with absolute values of coefficients
%
%   r = abs(p);
%
% r_i = abs(p_i), i.e. result is interval polynomial for interval polynomial input
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = abs(p.c);
