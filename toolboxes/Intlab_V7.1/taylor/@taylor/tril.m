function a = tril(a,k)
%TRIL         Implements  tril(a,k)  for Taylor
%
%   c = tril(a,k)
%
% functionality as Matlab function tril for matrices
%

% written  05/21/09     S.M. Rump
%

  if nargin==1
    k = 0;
  end

  index = reshape(1:prod(a.size),a.size);
  index = find(tril(index,k)==0);
  a.t(:,index) = 0;
