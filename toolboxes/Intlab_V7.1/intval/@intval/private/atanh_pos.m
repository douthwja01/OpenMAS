function y = atanh_pos(x,rnd)
% rigorous atanh for nonnegative double vector x<1 rounded corresponding to rnd
% for internal use in rigorous atanh, acoth

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  major revision
%

  y = x;

  % atanh(x) = log( 1 + 2*x/(1-x) ) / 2

  index = ( x<0.33 );                    % 0 <= x < .33
  if any(index)
    X = x(index);
    setround(-rnd)
    E = 1 - X;
    setround(rnd)
    E = 2*X./E;
    y(index) = log_1( E , rnd ) / 2;     % 0 <= E < 1
  end

  index = ~index;
  if any(index)                          % .33 <= x <= 1
    X = x(index);
    setround(-rnd)
    E = 1 - X;
    setround(rnd)
    E = 1 + 2*X./E;
    y(index) = log_rnd( E , rnd ) / 2;
  end
