function res = ne(a,b)
%NE           Implements  a ~= b  for hessians, compares only a.x and b.x
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/21/09     S.M. Rump  typo
%

  if ~isa(a,'hessian')
    res = ( a~=b.x );
  elseif ~isa(b,'hessian')
    res = ( a.x~=b );
  else
    res = ( a.x~=b.x );
  end
