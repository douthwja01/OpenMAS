function C = power(A,n)
%POWER        Implements  A .^ n  for long numbers  (same as A^n)
%
% Input A may long number or (column) vector, n is integer
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  C = A^n;
