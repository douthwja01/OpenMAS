function r = issparse(c)
%ISSPARSE     Returns 1 if c is sparse
%
%  r = issparse(c);
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if c.complex
    r = issparse(c.mid);
  else
    r = issparse(c.inf);
  end
