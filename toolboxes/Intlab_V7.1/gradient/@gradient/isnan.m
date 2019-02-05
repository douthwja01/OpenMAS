function r = isnan(a)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    Matlab flaw for empty sparse corrected
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = reshape( isnan(a.x(:)) | any(isnan(a.dx),2) , size(a.x) );
