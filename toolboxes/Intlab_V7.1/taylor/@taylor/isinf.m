function r = isinf(a)
%ISINF        Array of 1's for inf components
%
%   r = isinf(a)
%

% written  05/21/09     S.M. Rump
%

  r = reshape(any(isinf(a.t)),a.size);
