function r = isinf(a)
%ISINF        Array of 1's for inf components
%
%   r = isinf(a)
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    Matlab flaw for empty sparse corrected
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( isinf(a.x(:)) | any(isinf(a.dx),2) , size(a.x) );
