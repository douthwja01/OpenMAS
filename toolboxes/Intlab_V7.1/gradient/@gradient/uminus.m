function a = uminus(a)
%UMINUS       Gradient unary minus  - a
%

% written  10/16/98     S.M. Rump
% modified 12/18/02     S.M. Rump  performance improved
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = -a.x;
  a.dx = -a.dx;
