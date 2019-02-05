function res = all(a,dim)
%ANY          Like Matlab function "all" for Taylor
%
%Call
%
%   L = all(A)
%   L = all(A,dim)
%
%Same functionality as Matlab/all for intval quantity A
%

% written  05/21/09     S.M. Rump
%

  if nargin==1
    res = any(reshape(all(a.t),a.size));
  else
    res = any(reshape(all(a.t),a.size),dim);
  end
