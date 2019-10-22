function a = mid(a)
%MID          Taylor midpoint
%
%  r = mid(a)
%

% written  05/21/09     S.M. Rump
%

  if isa(a.t,'intval')
    a.t = mid(a.t);
  end
