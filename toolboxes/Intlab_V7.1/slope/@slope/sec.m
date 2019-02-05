function a = sec(a)
%SEC          Slope secant sec(a)
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = 1 ./ cos(a);
