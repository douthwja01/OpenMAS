function A = abs(A)
%ABS          Absolute value of long (refers only to midpoint)
%
%  C = abs(A)
%
%For computation with error term and  A = midA +/- err, on output
%  C = |midA| +/- err .
%

% written  11/06/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 07/25/05     S.M. Rump  improved performance
%

  A.sign = abs(A.sign);
