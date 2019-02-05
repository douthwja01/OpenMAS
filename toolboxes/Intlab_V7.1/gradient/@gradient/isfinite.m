function r = isfinite(a)
%ISFINITE     Array of 1's for finite components
%
%   r = isfinite(a)
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = isfinite(a.x) & reshape(all(isfinite(a.dx),2),size(a.x));
