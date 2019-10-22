function c = diag(a,k)
%DIAG         Implements  diag(a,k)  for intervals
%
%   c = diag(a,k)
%
% functionality as Matlab function diag for matrices
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
    c.mid = diag(a.mid,k);
    c.rad = diag(a.rad,k);
  else
    c.inf = diag(a.inf,k);
    c.sup = diag(a.sup,k);
  end
