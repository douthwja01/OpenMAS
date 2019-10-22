function x = sqr(x)
%SQR          Elementwise square of double (for completeness for intval, gradient, slope)
%
%  res = sqr(x)
%

% written  12/06/95     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  x = x.*x;
