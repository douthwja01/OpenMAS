function res = le(a,b)
%LE           Implements  a <= b  for Taylor, compares only a.x and b.x
%

% written  05/21/09     S.M. Rump
%

  if ~isa(a,'taylor')
    res = ( a<=reshape(b.t(1,:),b.size) );
  elseif ~isa(b,'taylor')
    res = ( reshape(a.t(1,:),a.size)<=b );
  else
    res = ( reshape(a.t(1,:),a.size)<=reshape(b.t(1,:),b.size) );
  end
