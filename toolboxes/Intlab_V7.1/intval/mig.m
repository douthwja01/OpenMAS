function res = mig(a)
%MIG          Mignitude of point matrix (for completeness)
%
%   res = mig(a)
%

% written  11/23/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  res = abs(a);
  