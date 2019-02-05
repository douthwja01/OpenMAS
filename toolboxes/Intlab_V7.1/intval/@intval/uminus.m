function a = uminus(a)
%UMINUS       Implements  -a  for intervals
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged, improved performance
%

  if a.complex
    a.mid = -a.mid;
  else
    asup = a.sup;
    a.sup = -a.inf;
    a.inf = -asup;
  end
