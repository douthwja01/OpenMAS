function a = full(a)
%FULL         Convert gradient to full
%

% written  03/06/04     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = full(a.x);
  a.dx = full(a.dx);
