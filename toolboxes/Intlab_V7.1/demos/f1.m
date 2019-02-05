function y = f1(x)
  y = x;
  y(1) = .5*sin(x(1)*x(2)) - x(2)/(4*pi) - x(1)/2;
  y(2) = (1-1/(4*pi))*(exp(2*x(1))-exp(1)) + exp(1)*x(2)/pi - 2*exp(1)*x(1);