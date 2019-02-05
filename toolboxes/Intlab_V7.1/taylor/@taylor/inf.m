function a = inf(a)
%INF          Taylor infimum
%
%  r = inf(a)
%

% written  05/21/09     S.M. Rump
%

  if isa(a.t,'intval')
    a.t = inf(a.t);
  end
