function r = isnan(a)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(a)
%

% written  09/29/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = isnan(a.sign) | any(isnan(a.mantissa),2) | isnan(a.exponent);
