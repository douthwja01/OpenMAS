function a = inf(a)
%INF          Hessian infimum
%
%  r = inf(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isa(a.x,'intval')
    a.x = inf(a.x);
    a.dx = inf(a.dx);
    a.hx = inf(a.hx);
  end
