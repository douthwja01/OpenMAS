function a = inf(a)
%INF          Gradient infimum
%
%  r = inf(a)
%

% written  10/16/98     S.M. Rump
% modified 12/18/02     S.M. Rump  improved performance
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isa(a.x,'intval')
    a.x = inf(a.x);
    a.dx = inf(a.dx);
  end
