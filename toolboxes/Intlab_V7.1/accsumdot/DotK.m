function res = DotK(x,y,K)
%DOTK         Dot product 'as if' computed in K-fold (quadruple) precision
%
%   res = DotK(x,y,K)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in K-fold precision. For K=inf, AccSum will be used for summation.
%  Note that the result need to be faithful if underflow occurs.
%
%Implements algorithm DotK from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires (12K+1)n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input, K=inf
%

  if ~isreal(x)
    if ~isreal(y)                       % both x and y complex
      x = x(:);
      y = y(:);
      res = complex( DotK([real(x);imag(x)],[real(y);-imag(y)],K) , ...
                     DotK([real(x);imag(x)],[imag(y);real(y)],K) );
    else                                % x complex, y real
      res = complex( DotK(real(x),y,K) , DotK(imag(x),y,K) );
    end
    return
  else                                  % x real
    if ~isreal(y)                       % y complex
      res = complex( DotK(x,real(y),K) , DotK(x,imag(y),K) );
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

  n = length(x);
  r = zeros(2*n,1);
  [p,r(1)] = TwoProduct(x(1),y(1));
  for i=2:n
    [h,r(i)] = TwoProduct(x(i),y(i));
    [p,r(n+i-1)] = TwoSum(p,h);
  end
  r(2*n) = p;
  if K==inf
    res = AccSum(r);
  else
    res = SumK(r,K-1);
  end
  
  if rndold
    setround(rndold)
  end
  