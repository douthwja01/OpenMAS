function a = imag(a)
%IMAG         imaginary part of hessians
%
%   c = imag(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = imag(a.x);
  a.dx = imag(a.dx);
  a.hx = imag(a.hx);
