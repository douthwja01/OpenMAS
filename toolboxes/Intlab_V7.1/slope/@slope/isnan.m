function r = isnan(u)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(u)
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( any(isnan(u.r),2) | any(isnan(u.s),2) , u.size );
