function res = Sum2(p)
%SUM2         Summation 'as if' computed in 2-fold (quadruple) precision
%
%   res = Sum2(p)
%
%On return, res approximates sum(p) with accuracy as if computed 
%  in 2-fold precision. Input vector p may be single or double precision.
%
%Implements algorithm Sum2 from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires 7n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

% Implementation in the paper
%   for i=2:length(p)
%     [p(i),p(i-1)] = TwoSum(p(i),p(i-1));
%   end
%   res = sum(p(1:end-1)) + p(end);
  

  if ~isreal(p)
    res = complex(Sum2(real(p)),Sum2(imag(p)));
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  % Implementation without overwriting p and expanded TwoSum
  e = 0;
  s = p(1);
  for i=2:length(p)
    x = p(i) + s;                   % [x,y] = TwoSum(p_i,s)
    z = x - p(i);
    y = ( p(i) - (x-z) ) + (s-z);
    e = e + y;
    s = x;
  end
  res = s + e;
  
  if rndold
    setround(rndold)
  end
  