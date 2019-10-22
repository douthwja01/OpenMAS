function a = transpose(a)
%TRANSPOSE    Implements  a.'  for Hessians
%
%   c = a.'
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  [m n] = size(a.x);
  a.x = a.x.';
  if m*n~=1
    index = reshape( 1:(m*n) , m , n )';
    a.dx = a.dx( : , index(:) );
    a.hx = a.hx( : , index(:) );
  end
