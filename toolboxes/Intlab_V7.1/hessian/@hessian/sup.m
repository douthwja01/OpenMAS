function a = sup(a)
%SUP          Hessian supremum
%
%  r = sup(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isa(a.x,'intval')
    a.x = sup(a.x);
    a.dx = sup(a.dx);
    a.hx = sup(a.hx);
  end
