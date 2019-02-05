function r = issparse(c)
%ISSPARSE     Returns 1 if c is sparse
%
%  r = issparse(c)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    functionality changed: true if .dx is sparse
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = issparse(c.dx);
