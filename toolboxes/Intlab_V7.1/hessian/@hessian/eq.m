function res = eq(a,b)
%EQ           Implements  a == b  for hessians, compares only a.x and b.x
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isa(a,'hessian')
    res = ( a==b.x );
  elseif ~isa(b,'hessian')
    res = ( a.x==b );
  else
    res = ( a.x==b.x );
  end
