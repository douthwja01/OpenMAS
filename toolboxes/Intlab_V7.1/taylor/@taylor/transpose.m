function a = transpose(a)
%TRANSPOSE    Implements  a.'  for Taylor
%
%   c = a.'
%

% written  05/21/09     S.M. Rump
%

  index = reshape(1:prod(a.size),a.size)';
  a.size = fliplr(a.size);
  a.t = a.t(:,index);
