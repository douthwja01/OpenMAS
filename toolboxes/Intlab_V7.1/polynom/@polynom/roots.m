function r = roots(p)
%ROOTS        Approximations to roots of univariate non-interval polynomial
%
%   r = roots(p)
%

% written  09/14/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if size(p.e,2)>1
    error('roots only for univariate polynomials')
  end
  if isa(p.c,'intval')
    error('roots only for non-interval polynomials')
  end
  r = roots(p.c);
