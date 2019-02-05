function y = testfuntaylor(x)
%TESTFUNTAYLOR   A test function for the Taylor toolbox
%

% written  05/29/09     S.M. Rump
% modified 08/22/12     S.M. Rump  Constant pi
%

  if isintval(x)
    % Pi = 4*atan(intval(1));
    Pi = intval('pi');
  else
    Pi = pi;
  end
  y = sin(Pi*x)-sin(x);
