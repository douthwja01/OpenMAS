function c = mldivide(a,b)
%MLDIVIDE     Implements  a \ b
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  c = verifylss(a,b);
