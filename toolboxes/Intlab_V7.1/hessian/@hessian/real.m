function a = real(a)
%REAL         real part of hessians
%
%   c = real(a)
%
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = real(a.x);
  a.dx = real(a.dx);
  a.hx = real(a.hx);
