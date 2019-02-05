function r = mpower(a,b)
%MPOWER       Implements  a ^ b  for gradients
%
% either a and b are (interval) scalar or, b is integer
%

% written  10/16/98     S.M. Rump
% modified 12/30/98     S.M. Rump  improved performance and even exponent
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  % for scalars use .^ with improved diameter for even exponent
  if ( prod(size(a))==1 ) & ( prod(size(b))==1 )
    r = a .^ b ;
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isa(b,'double') & isreal(b) & prod(size(b))==1 & b==round(b)
    [m n] = size(a);
    if m~=n
      error('gradient mpower of non-square matrix')
    end
    if b==0
        r = typeadj( gradient(eye(size(a))) , typeof(a) );
    elseif b<0                  % b is negative integer
      error('negative gradient matrix power')
    else                        % b is integer
      b = b - 1;                % b is at least 1
      r = a;
      while b>0
        if mod(b,2)==1
          r = r*a;
        end
        b = floor(b/2);
        if b~=0
          a = a*a;
        end
      end
    end
  else
    error('invalid call of gradient mpower ^')
  end
  
  if rndold
    setround(rndold)
  end
