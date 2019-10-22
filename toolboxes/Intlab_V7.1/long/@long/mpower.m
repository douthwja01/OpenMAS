function C = mpower(A,n)
%MPOWER       Implements  A ^ n  for long numbers
%
% Input A may long number or vector, n is integer
%

% written  12/30/98     S.M. Rump
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

  if isa(n,'double') & isreal(n) & prod(size(n))==1 & n==round(n)
    if n==0
      C = long(1);
    else                        % n is integer
      n_is_negative = ( n<0 );
      n = abs(n) - 1;           % abs(n) is at least 1
      C = A;
      while n>0
        if mod(n,2)==1
          C = C*A;
        end
        n = floor(n/2);
        if n~=0
          A = A*A;
        end
      end
      if n_is_negative
        C = 1/C;
      end
    end
  else
    error('invalid call of mpower ^')
  end
  
  if rndold
    setround(rndold)
  end
