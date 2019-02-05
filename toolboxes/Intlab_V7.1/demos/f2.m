function  y = f2(x)
  y = x;
  cPi = typeadj( midrad(3.141592653589793,1e-15) , typeof(x) );
  y(1) = .5*sin(x(1)*x(2)) - x(2)/(4*cPi) - x(1)/2;
  y(2) = (1-1/(4*cPi))*(exp(2*x(1))-exp(1)) + exp(1)*x(2)/cPi - 2*exp(1)*x(1);