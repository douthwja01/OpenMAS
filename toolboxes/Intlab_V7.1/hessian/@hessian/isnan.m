function r = isnan(a)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( isnan(a.x(:)') | any(isnan(a.dx),1) | any(isnan(a.hx),1) , size(a.x) );
