function c = conj(a)
%CONJ         Implements complex conjugate  conj(a)  for intervals
%
%   c = conj(a)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  c = a;
  if a.complex
    c.mid = conj(a.mid);
  end
