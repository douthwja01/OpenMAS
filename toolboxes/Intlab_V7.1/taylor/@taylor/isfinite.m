function r = isfinite(a)
%ISFINITE     Array of 1's for finite components
%
%   r = isfinite(a)
%

% written  05/21/09     S.M. Rump
%

  r = reshape(all(isfinite(a.t)),a.size);
