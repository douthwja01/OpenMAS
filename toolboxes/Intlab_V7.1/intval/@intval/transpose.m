function a = transpose(a)
%TRANSPOSE    Implements  a.'  for intervals (not conjugate!)
%
%  c = a.'
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged, improved performance
%

  if a.complex
    a.mid = a.mid.';
    a.rad = a.rad.';
  else
    a.inf = a.inf.';
    a.sup = a.sup.';
  end
