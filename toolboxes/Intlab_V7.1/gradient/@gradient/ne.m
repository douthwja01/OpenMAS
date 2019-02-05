function res = ne(a,b)
%NE           Implements  a ~= b  for gradients, compares only a.x and b.x
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/21/09     S.M. Rump  typo
%

  if ~isa(a,'gradient')
    res = ( a~=b.x );
  elseif ~isa(b,'gradient')
    res = ( a.x~=b );
  else
    res = ( a.x~=b.x );
  end
