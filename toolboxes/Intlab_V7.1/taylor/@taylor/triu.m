function a = triu(a,k)
%TRIU         Implements  triu(a,k)  for Taylor
%
%   c = triu(a,k)
%
% functionality as Matlab function triu for matrices
%

% written  05/21/09     S.M. Rump
%

  if nargin==1
    k = 0;
  end

  index = reshape(1:prod(a.size),a.size);
  index = find(triu(index,k)==0);
  a.t(:,index) = 0;
