function Y = full(X)
%FULL         type cast to full interval matrix
%
%   Y = full(X)
%

% written  06/21/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  Y = X;

  if X.complex
    Y.mid = full(X.mid);
    Y.rad = full(X.rad);
  else
    Y.inf = full(X.inf);
    Y.sup = full(X.sup);
  end
