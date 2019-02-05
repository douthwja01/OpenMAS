function r = isfinite(u)
%ISFINITE     Array of 1's for finite components
%
%   r = isfinite(u)
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( all(isfinite(u.r),2) & all(isfinite(u.s),2) , u.size );
