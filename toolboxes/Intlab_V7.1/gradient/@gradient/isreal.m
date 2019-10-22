function r = isreal(a)
%ISREAL       returns 1 if c is real
%
%   r = isreal(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of derivative
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = isreal(a.x) & isreal(a.dx);
