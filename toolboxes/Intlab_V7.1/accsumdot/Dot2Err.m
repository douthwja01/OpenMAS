function [res,err] = Dot2Err(x,y)
%DOT2Err      Dot product 'as if' computed in 2-fold (quadruple) precision with error bound
%
%   [res,err] = Dot2Err(x,y)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in 2-fold precision, and (also valid in case of underflow):
%
%  | x'*y - res | <= err .
%
%Implements algorithm Dot2Err from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires 27n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ( ~isreal(x) ) | ( ~isreal(y) )
    error('Dot2Err for real input only')
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  n = length(x);
  if isa(x,'double'), prec='double'; else prec='single'; end
  Eps = eps(prec)/2;
  eta = succ(0*eps);
  if 2*n*eps>=1
    error('Dot2Err: vector length too large')
  end
  [p,s] = TwoProduct(x(1),y(1));
  e = abs(s);
  for i=2:length(x)
    [h,r] = TwoProduct(x(i),y(i));
    [p,q] = TwoSum(p,h);
    t = q + r;
    s = s + t;
    e = e + abs(t);
  end
  res = p + s;
  delta = n*Eps/(1-2*n*Eps);
  alpha = Eps*abs(res) + (delta*e + 3*eta/Eps);
  err = alpha / (1-2*Eps);
  
  if rndold
    setround(rndold)
  end
  