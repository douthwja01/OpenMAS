function res = qdist(a,b)
%QDIST        Implements  q(a,b)  metrical distance
%  Name  qdist  is used to avoid ambiguities with variable  q
%  This functions for non-interval input only for completeness
%
%     res = qdist(a,b)
%
% for real input          abs(a-b)
% for complex input       qdist(real(a),real(b)) + qdist(imag(a),imag(b))
%

% written  03/23/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isreal(a) & isreal(b)
    res = abs(a-b);
  else
    res = abs(real(a)-real(b)) + abs(imag(a)-imag(b));
  end
