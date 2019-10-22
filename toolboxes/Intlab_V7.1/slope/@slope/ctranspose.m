function u = ctranspose(a)
%TRANSPOSE    Implements  a'  for slopes
%
%   u = a'
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if length(a.size)>2
    error('transpose not defined for multi-dimensional arrays')
  end

  u = a;

  u.size = fliplr(u.size);
  index = reshape( 1:prod(a.size) , a.size )';
  index = index(:);
  u.r = u.r( index , : );
  u.s = u.s( index , : );
