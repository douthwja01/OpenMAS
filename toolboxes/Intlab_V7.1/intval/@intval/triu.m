function c = triu(a,k)
%TRIU         Implements  triu(a,k)  for intervals
%
%   c = triu(a,k)
%
% functionality as Matlab function triu for matrices
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
    c.mid = triu(a.mid,k);
    c.rad = triu(a.rad,k);
  else
    c.inf = triu(a.inf,k);
    c.sup = triu(a.sup,k);
  end
