function a = full(a)
%FULL         Convert hessian first and second derivative to full
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = full(a.x);
  a.dx = full(a.dx);
  a.hx = full(a.hx);