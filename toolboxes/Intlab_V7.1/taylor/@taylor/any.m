function res = any(a,dim)
%ANY          Like Matlab function "any" for Taylor
%
%Call
%
%   L = any(A)
%   L = any(A,dim)
%
%Same functionality as Matlab/any for intval quantity A
%

% written  05/21/09     S.M. Rump
%

  if nargin==1
    res = any(reshape(any(a.t),a.size));
  else
    res = any(reshape(any(a.t),a.size),dim);
  end
