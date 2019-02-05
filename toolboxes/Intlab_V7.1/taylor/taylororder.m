function k = taylororder
%TAYLORORDER  Order of Taylor evaluations
%
%   k = taylororder
%
%In the current initialization of the Taylor package, Taylor coefficients
%up to order k will be calculated. Result is -1 if Taylor package is not
%yet initialized.
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_TAYLOR_ORDER = getappdata(0,'INTLAB_TAYLOR_ORDER');
  if isempty(INTLAB_TAYLOR_ORDER)   % Taylor package not yet initialized
    k = -1;
  else
    k = INTLAB_TAYLOR_ORDER;
  end
