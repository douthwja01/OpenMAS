function a = uminus(a)
%UMINUS       Hessian unary minus  - a
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = -a.x;
  a.dx = -a.dx;
  a.hx = -a.hx;
