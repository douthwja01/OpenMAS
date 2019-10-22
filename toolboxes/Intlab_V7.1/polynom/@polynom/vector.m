function r = vector(p);
%VECTOR       (row) vector of coefficients of univariate polynomial
%
%   r = vector(p);
%

% written  08/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if size(p.e,2)>1
    error('input multivariate polynomial')
  else
    r = p.c;
  end
