function c = trace(a)
%TRACE        Implements  trace(a)  for Taylor
%
%   c = trace(a)
%

% written  05/21/09     S.M. Rump
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  c = sum( diag( a ) );
  
  if rndold
    setround(rndold)
  end
