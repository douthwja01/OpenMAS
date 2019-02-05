function a = sup(a)
%SUP          Taylor supremum
%
%  r = sup(a)
%

% written  05/21/09     S.M. Rump
%

  if isa(a.t,'intval')
    a.t = sup(a.t);
  end
