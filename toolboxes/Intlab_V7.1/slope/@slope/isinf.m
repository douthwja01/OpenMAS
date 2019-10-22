function r = isinf(u)
%ISINF        Array of 1's for inf components
%
%   r = isinf(u)
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( any(isinf(u.r),2) | any(isinf(u.s),2) , u.size );
