function res = ufp(x)
%UFP          unit in the first place (ufp) of real (vector, matrix) x
%
%   res = ufp(x)
%
%

% written  10/18/08     S.M. Rump
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  [f,e] = log2(abs(x));                 % the easy way
  res = pow2(floor(2*f)/2,e);

  if rndold
    setround(rndold)
  end
