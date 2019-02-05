function C = times(A,B)
%TIMES        Long multiplication  A .* B  (same as A * B)
%
%Multiplication always entrywise
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  C = A * B;
