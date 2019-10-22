function r = power(a,b)
%MPOWER       Hessian elementwise power  a .^ b
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 05/14/09     S.M. Rump  NaN^0
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isa(b,'double') & isreal(b) & prod(size(b))==1 & b==round(b)
    if b==0
      r = typeadj( ones(size(a)) , typeof(a) );
      r(isnan(a)) = NaN;
    else                        % b is integer
      b_is_negative = b<0;
      % check b is even to ensure result is nonnegative
      if b==2*floor(b/2)
        b = b/2;
        b_is_even = 1;
      else
        b_is_even = 0;
      end
      b = abs(b) - 1;           % abs(b) is at least 1
      r = a;
      while b>0
        if mod(b,2)==1
          r = r.*a;
        end
        b = floor(b/2);
        if b~=0
          a = sqr(a);
        end
      end
      if b_is_even
        r = sqr(r);
      end
      if b_is_negative
        r = 1./r;
      end
    end
  else
    r = exp( b .* log(a) );
  end
  
  if rndold
    setround(rndold)
  end
