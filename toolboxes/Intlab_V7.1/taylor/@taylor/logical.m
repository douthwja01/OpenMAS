function res = logical(a)
%LOGICAL      Like Matlab function "logical" for Taylor
%
%Call
%
%   L = logical(A)
%
%Same functionality as Matlab/logical for intval quantity A
%

% written  05/21/09     S.M. Rump
%

  res = reshape(logical(any(a.t)),a.size);
