function c = band(a,p,q)
%BAND         Extract band from matrix a, lower bandwidth p, upper bandwidth q
%   if parameter q is omitted, q:=p
%
%   c = band(a,p,q)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin<3
    q = p;
  end

  c = a;

  if a.complex
    c.mid = tril(triu(a.mid,-p),q);
    c.rad = tril(triu(a.rad,-p),q);
  else
    c.inf = tril(triu(a.inf,-p),q);
    c.sup = tril(triu(a.sup,-p),q);
  end
