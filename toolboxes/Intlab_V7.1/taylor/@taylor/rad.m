function a = rad(a)
%RAD          Taylor radius
%
%  r = rad(a)
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  a.t = rad(a.t);
  if isequal(a.t,0)
    K1 = getappdata(0,'INTLAB_TAYLOR_ORDER') + 1;
    a.t = zeros(K1,prod(a.size));
  end
