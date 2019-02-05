function d = discriminant(p)
%DISCRIMINANT Discriminant of a univariate polynomial p
%
%   d = discriminant(p)
%
%

% written  02/04/04     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  d = det(sylvester(p,p'));
