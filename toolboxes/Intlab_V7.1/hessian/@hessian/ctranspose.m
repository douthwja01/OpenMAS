function a = ctranspose(a)
%CTRANSPOSE   Implements  a'  for Hessians
%
%   c = a'
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isreal(a.x)
    error('Hessian function '' only for real arguments')
  end

  [m n] = size(a.x);
  a.x = a.x';
  if m*n~=1
    index = reshape( 1:(m*n) , m , n )';
    a.dx = a.dx( : , index(:) );
    a.hx = a.hx( : , index(:) );
  end
