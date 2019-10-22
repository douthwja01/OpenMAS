function r = isnan(a)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(a)
%

% written  05/21/09     S.M. Rump
%

  r = reshape(any(isnan(a.t)),a.size);
