function a = abss(a)
%ABSS         Absolute value for non-intervals, same as abs (for completeness)
%
%  c = abss(a);
%
%Obsolete, replaced by mag (thanks to Arnold for proposing better naming).
%

% written  07/08/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = mag(a);
