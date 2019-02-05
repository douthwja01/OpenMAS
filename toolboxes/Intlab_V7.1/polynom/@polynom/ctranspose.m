function r = ctranspose(p)
%CTRANSPOSE   Implements  p'  for univariate polynomials (derivative)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = pderiv(p);
