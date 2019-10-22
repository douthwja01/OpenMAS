function r = isintval(a)
%ISINTVAL     Returns 1 if  a  is intval slope
%
%   r = isintval(a)
%
%For completeness:  returns 1 anyway because slopes components are always intval
%

% written  09/02/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = logical(1);
