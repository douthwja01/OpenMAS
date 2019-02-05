function r = isreal(c)
%ISREAL       returns 1 if c is real
%
%  r = isreal(c)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = ~c.complex;
