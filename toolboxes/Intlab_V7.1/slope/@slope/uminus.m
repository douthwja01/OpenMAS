function a = uminus(a)
%UMINUS       Slope unary minus  - a
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.r = - a.r;
  a.s = - a.s;
