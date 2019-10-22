function res = qdist(a,b)
%QDIST        Implements  q(a,b)  metrical distance
%  Name  qdist  is used to avoid ambiguities with variable  q
%  This functions for non-interval input only for completeness
%
%     res = qdist(a,b)
%
% for real input          abs(a.c-b.c)
% for complex input       qdist(real(a.c),real(b.c)) + qdist(imag(a.c),imag(b.c))
%

% written  10/04/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  a = polynom(a);
  b = polynom(b);

  if isreal(a.c) & isreal(b.c)
    res = abs(a.c-b.c);
  else
    res = abs(real(a.c)-real(b.c)) + abs(imag(a.c)-imag(b.c));
  end
  
  if rndold
    setround(rndold)
  end
