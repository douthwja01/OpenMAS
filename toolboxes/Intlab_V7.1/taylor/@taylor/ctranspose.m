function a = ctranspose(a)
%CTRANSPOSE   Implements  a'  for Taylor
%
%   c = a'
%

% written  05/21/09     S.M. Rump
%

  if ~isreal(a.t)
    error('taylor function '' only for real arguments')
  end

  index = reshape(1:prod(a.size),a.size)';
  a.size = fliplr(a.size);
  a.t = conj(a.t(:,index));
