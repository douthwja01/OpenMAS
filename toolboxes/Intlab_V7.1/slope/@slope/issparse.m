function r = issparse(u)
%ISSPARSE     Returns 1 if u is sparse
%
%  r = issparse(u)
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = issparse(u.r);
