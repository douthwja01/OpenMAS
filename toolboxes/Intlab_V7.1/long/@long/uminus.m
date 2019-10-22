function A = uminus(A)
%UMINUS       Long unary minus  -B
%
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  A.sign = -A.sign;
