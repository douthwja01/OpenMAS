function  y = demotest(x,radius)
%DEMOTEST     Test function for demoslope
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  y = x;
  if nargin==1
    radius = 1e-15;
  end
  cPi = typeadj( midrad(3.141592653589793,radius) , typeof(x) );
  y(1) = .5*sin(x(1)*x(2)) - x(2)/(4*cPi) - x(1)/2;
  y(2) = (1-1/(4*cPi))*(exp(2*x(1))-exp(1)) + exp(1)*x(2)/cPi - 2*exp(1)*x(1);
