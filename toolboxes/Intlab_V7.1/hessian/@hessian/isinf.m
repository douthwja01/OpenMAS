function r = isinf(a)
%ISINF        Array of 1's for inf components
%
%   r = isinf(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( isinf(a.x(:)') | any(isinf(a.dx),1) | any(isinf(a.hx),1) , size(a.x) );
