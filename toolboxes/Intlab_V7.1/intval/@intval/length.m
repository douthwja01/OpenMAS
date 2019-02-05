function m = length(a)
%LENGTH       Implements  length(a)  for intervals
%
%   m = length(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if a.complex
    m = length(a.mid);
  else
    m = length(a.inf);
  end
