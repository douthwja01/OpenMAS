function c = tril(a,k)
%TRIL         Implements  tril(a,k)  for intervals
%
%   c = tril(a,k)
%
% functionality as Matlab function tril for matrices
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin==1
    k = 0;
  end

  c = a;
  if a.complex
    c.mid = tril(a.mid,k);
    c.rad = tril(a.rad,k);
  else
    c.inf = tril(a.inf,k);
    c.sup = tril(a.sup,k);
  end
  