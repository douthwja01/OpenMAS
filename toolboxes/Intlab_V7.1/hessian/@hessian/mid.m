function a = mid(a)
%MID          Hessian midpoint
%
%  r = mid(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isa(a.x,'intval')
    a.x = mid(a.x);
    a.dx = mid(a.dx);
    a.hx = mid(a.hx);
  end
