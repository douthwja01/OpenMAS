function C = rdivide(A,B)
%RDIVIDE      Long right division  A ./ B  (same as A / B)
%
%Division always entrywise
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  C = A / B;
