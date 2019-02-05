function C = exp(X,p)
%EXP          Long exponential function
%
%  C = exp(X,p)
%
%for long number X. Input parameter p is optional; if specified,
%approximate accuracy of C approximately p decimals, otherwise that of X.
%Input X must be less than beta.
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  
  % convert decimal precision in beta-digits
  if nargin==2
    longprecision(p);
    p = ceil( p/log10(INTLAB_LONG_BETA) );
  else
    p = size(X.mantissa,2) + 1;
  end

  if any( X.exponent>1 )
    error('exponent too big for long exponential')
  end

  index = ( X.exponent>0 );
  if any( index )
    large = 1;
    E = zeros(size(X.sign));
    while any(index)
      E(index) = E(index) + 1;
      X(index) = X(index)/2;
      index = ( X.exponent>0 );
    end
  else
    large = 0;
  end

  C = 1 + X;
  T = X;
  i = 1;
  while 1
    i = i+1;
    T = T * X / i;
    if all( T.exponent<-p )
      if isequal( longinit('ErrorTerm',0) , 'WithErrorTerm' )
        C.error = errorupdate( 1 , C.error , 0 , 1 , 1 , -p );
      end
      break
    end
    C = C + T;
  end

  if large
    index = ( E~=0 );
    while any(index)
      C(index) = C(index)*C(index);
      E(index) = E(index) - 1;
      index = ( E~=0 );
    end
  end
  
  if rndold
    setround(rndold)
  end
