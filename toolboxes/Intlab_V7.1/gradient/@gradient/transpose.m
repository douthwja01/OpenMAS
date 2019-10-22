function a = transpose(a)
%TRANSPOSE    Implements  a.'  for gradients
%
%   c = a.'
%

% written  10/16/98     S.M. Rump
% modified 11/03/03     S.M. Rump  improved performance
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  [m n] = size(a.x);
  a.x = a.x.';
  if m*n~=1
    index = reshape( 1:(m*n) , m , n )';
    a.dx = a.dx( index , : );
  end
