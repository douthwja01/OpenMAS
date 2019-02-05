function r = isfinite(a)
%ISFINITE     Array of 1's for finite components
%
%   r = isfinite(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( isfinite(a.x(:)') & all(isfinite(a.dx),1) & all(isfinite(a.hx),1) , size(a.x) );
