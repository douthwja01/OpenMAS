function a = real(a)
%REAL         Real part of gradients
%
%   c = real(a)
%
%

% written  10/16/98     S.M. Rump
% modified 12/18/02     S.M. Rump  improved performance
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = real(a.x);
  a.dx = real(a.dx);
