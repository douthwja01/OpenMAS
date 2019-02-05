function c = spones(a);
%SPONES       Implements  spones(a)  for sparse interval matrix
%
%   c = spones(a)
%

% written  10/16/98     S.M. Rump
% modified 08/09/02     S.M. Rump  zero components inf/sup or mid/rad handled
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if a.complex
    c = spones(spones(a.mid)+spones(a.rad));
  else
    c = spones(spones(a.inf)+spones(a.sup));
  end
