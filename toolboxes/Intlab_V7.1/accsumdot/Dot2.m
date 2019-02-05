function res = Dot2(x,y)
%DOT2         Dot product 'as if' computed in 2-fold (quadruple) precision
%
%   res = Dot2(x,y)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in 2-fold precision.
%
%Implements algorithm Dot2 from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires 25n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ~isreal(x)
    if ~isreal(y)                       % both x and y complex
      x = x(:);
      y = y(:);
      res = complex( Dot2([real(x);imag(x)],[real(y);-imag(y)]) , ...
                     Dot2([real(x);imag(x)],[imag(y);real(y)]) );
    else                                % x complex, y real
      res = complex( Dot2(real(x),y) , Dot2(imag(x),y) );
    end
    return
  else                                  % x real
    if ~isreal(y)                       % y complex
      res = complex( Dot2(x,real(y)) , Dot2(x,imag(y)) );
      return
    end
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  [p,s] = TwoProduct(x(1),y(1));
  for i=2:length(x)
    [h,r] = TwoProduct(x(i),y(i));
    [p,q] = TwoSum(p,h);
    s = s + ( q + r );
  end
  res = p + s;
  
  if rndold
    setround(rndold)
  end
  