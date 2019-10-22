function r = mpower(a,b)
%MPOWER       Hessian power  a ^ n
%

% written  04/04/04     S.M. Rump
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
    if b<0
      error('negative exponent in mpower')
    end
    [m n] = size(a);
    if m~=n
      error('hessian mpower of non-square matrix')
    end
    if b==0
      if issparse(a)
        r = typeadj( hessian(speye(size(a))) , typeof(a) );
      else
        r = typeadj( hessian(eye(size(a))) , typeof(a) );
      end
    else                        % b is integer, at least 1
      b = b - 1;
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
    error('invalid call of hessian mpower ^')
  end
  
  if rndold
    setround(rndold)
  end
