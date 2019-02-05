function p = uminus(p)
%UMINUS       Polynomial unary minus  - p
%

% written  08/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged, improved performance
%

  p.c = -p.c;
