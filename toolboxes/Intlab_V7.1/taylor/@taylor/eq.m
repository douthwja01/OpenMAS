function res = eq(a,b)
%EQ           Implements  a == b  for Taylor, compares only function value, not derivatives
%

% written  05/21/09     S.M. Rump
%

  if ~isa(a,'taylor')
    res = ( a==reshape(b.t(1,:),b.size) );
  elseif ~isa(b,'taylor')
    res = ( reshape(a.t(1,:),a.size)==b );
  else
    res = ( reshape(a.t(1,:),a.size)==reshape(b.t(1,:),b.size) );
  end
