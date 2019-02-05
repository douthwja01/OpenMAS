function p = intval(p);
%INTVAL       Type conversion polynom to intval polynom
%
%   r = intval(p);
%

% written  08/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = intval(p.c);
