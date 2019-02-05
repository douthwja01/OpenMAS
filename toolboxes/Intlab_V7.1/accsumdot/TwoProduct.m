function [x,y] = TwoProduct(a,b)
%TWOPRODUCT   Error-free transformation of a+b into x*y with x=fl(a*b)
%
%   [x,y] = TwoProduct(a,b)
%
%On return, x+y=a*b and x=fl(a*b) provided no over- or underflow occurs .
%Input a,b may be vectors or matrices as well, in single or double precision.
%
%Follows G.W. Veltkamp, see T.J. Dekker: A floating-point technique for 
%  extending the available precision, Numerische Mathematik 18:224-242, 1971.
%Requires 17 flops for scalar input.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, underflow
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  x = a.*b;
  if any(~isfinite(x(:))) 
    error('overflow occurred in TwoProduct')
  end
  if isa(x,'double'), alpha=realmin('double'); else alpha=realmin('single'); end
  absx = abs(x(:));
  if any( ( absx<alpha ) & ( absx~=0 ))
    error('underflow occurred in TwoProduct')
  end
  [ah,al] = Split(a);
  [bh,bl] = Split(b);
  y = al.*bl - ( ( ( x - ah.*bh ) - al.*bh ) - ah.*bl );
  
  if rndold
    setround(rndold)
  end
