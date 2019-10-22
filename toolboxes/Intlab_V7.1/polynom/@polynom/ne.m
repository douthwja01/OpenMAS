function res = ne(p,q)
%NE           Implements  p ~= q  for polynomials
%
%Result 1 iff p and q are mathematically not identical
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p = removevars(permvars(polynom(p),'lex'));
  q = removevars(permvars(polynom(q),'lex'));
  res = ~( isequal(p.e,q.e) & isequal(p.c,q.c) & isequal(p.v,q.v) );
