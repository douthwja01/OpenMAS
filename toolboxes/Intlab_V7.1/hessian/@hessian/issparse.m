function r = issparse(c)
%ISSPARSE     Returns 1 if c.dx or c.hx is sparse
%
%  r = issparse(c)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = issparse(c.dx) | issparse(c.hx);
