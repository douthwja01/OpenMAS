function c = band(a,p,q)
%BAND         Extract band from matrix a, lower bandwidth p, upper bandwidth q
%   if parameter q is omitted, q:=p
%
%   c = band(a,p,q)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin<3
    q = p;
  end

  c = tril(triu(a,-p),q);
