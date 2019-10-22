function C = plus(A,B)
%PLUS         Long addition  A + B
%
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
  INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  A = long(A);
  B = long(B);

  % check dimensions
  nA = size(A.mantissa,1);
  nB = size(B.mantissa,1);
  if ( nA==1 ) & ( nB~=1 )
    Ones = ones(nB,1);
    A.sign = Ones*A.sign;
    A.mantissa = Ones*A.mantissa;
    A.exponent = Ones*A.exponent;
    if INTLAB_LONG_ERROR
      A.error.mant = Ones*A.error.mant;
      A.error.exp = Ones*A.error.exp;
    end
    nA = nB;
  end
  if ( nA~=1 ) & ( nB==1 )
    Ones = ones(nA,1);
    B.sign = Ones*B.sign;
    B.mantissa = Ones*B.mantissa;
    B.exponent = Ones*B.exponent;
    if INTLAB_LONG_ERROR
      B.error.mant = Ones*B.error.mant;
      B.error.exp = Ones*B.error.exp;
    end
    nB = nA;
  end
  if ( nA~=nB )
    error('incompatible operands for long plus')
  end

  % set sign
  C.sign = A.sign;

  C.exponent = max( A.exponent , B.exponent );
  A = adjustexp(A,C.exponent-A.exponent);
  B = adjustexp(B,C.exponent-B.exponent);

  % indices of different sign ("minus-components")
  indexminus = ( A.sign~=B.sign );

  % add mantissas
  precA = size(A.mantissa,2);
  precB = size(B.mantissa,2);
  if precA>=precB
    C.mantissa = A.mantissa;
    C.mantissa(:,1:precB) = A.mantissa(:,1:precB) + B.mantissa;
    if any(indexminus)
      C.mantissa(indexminus,1:precB) = ...
         A.mantissa(indexminus,1:precB) - B.mantissa(indexminus,:);
    end
    if INTLAB_LONG_ERROR
      C.error = errorupdate( 1 , A.error , 0 , 1 , B.error , 0 );
    end
  else
    C.sign = B.sign;
    C.mantissa = B.mantissa;
    C.mantissa(:,1:precA) = B.mantissa(:,1:precA) + A.mantissa;
    if any(indexminus)
      C.mantissa(indexminus,1:precA) = ...
         B.mantissa(indexminus,1:precA) - A.mantissa(indexminus,:);
    end
    if INTLAB_LONG_ERROR
      C.error = errorupdate( 1 , A.error , 0 , 1 , B.error , 0 );
    end
  end

  % special treatment of "minus-components": A and B of different sign
  % prenormalize minus-components
  C.mantissa(indexminus,:) = C.mantissa(indexminus,:) + INTLAB_LONG_BETA - 1;
  C.mantissa(indexminus,end) = C.mantissa(indexminus,end) + 1;

  C = normalize(C);

  % adjust minus-components
  index = ( C.mantissa(:,1)>=INTLAB_LONG_BETA );
  indexpos = indexminus & index;
  indexneg = indexminus & ~index;

  % adjust minus-components w/o change of sign
  if any(indexpos)   % no change of sign, renormalize
    C.mantissa(indexpos,1) = C.mantissa(indexpos,1) - INTLAB_LONG_BETA;
  end

  % adjust minus-components with change of sign
  if any(indexneg)   % change of sign
    C.sign(indexneg) = -C.sign(indexneg);
    C.mantissa(indexneg,:) = INTLAB_LONG_BETA-1 - C.mantissa(indexneg,:);
    C.mantissa(indexneg,end) = C.mantissa(indexneg,end) + 1;

    %normalize again (plus 1 might cause last digit equal to beta)
    while any(any( C.mantissa(indexneg,:)>=INTLAB_LONG_BETA ))
      q = floor( C.mantissa(indexneg,:)/INTLAB_LONG_BETA );
      C.mantissa(indexneg,:) = C.mantissa(indexneg,:) - q*INTLAB_LONG_BETA;
      C.mantissa(indexneg,1:end-1) = C.mantissa(indexneg,1:end-1) + q(:,2:end);
    end
  end

  % Normalize first digit
  C = normalizefirst(C);

  % Omit leading zeros
  C = omitleadingzeros(C);
  if INTLAB_LONG_ERROR
    C.error = errornormalize(C.error);
  end

  % If necessary, cut length
  precC = size(C.mantissa,2);
  if precC>INTLAB_LONG_PRECISION
    if INTLAB_LONG_ERROR
      C.error = errorupdate( 1 , C.error , 0 , ...
                   1 , any( C.mantissa(:,INTLAB_LONG_PRECISION+1:end),2 ) , ...
                   C.exponent-INTLAB_LONG_PRECISION );
    end
    C.mantissa = C.mantissa(:,1:INTLAB_LONG_PRECISION);
  end
  if INTLAB_LONG_ERROR
    C.error = errornormalize(C.error);
  else
    C.error.mant = 0;
    C.error.exp = 0;
  end

  C = class(C,'long');
  
  if rndold
    setround(rndold)
  end


function A = adjustexp(A,d)
%adjust exponent of A by d beta-digits: add d leading (beta-) zeros
%subject to maximum precision (careful: d is vector)
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  dnonzero = ( d~=0 );
  if any(dnonzero)

    % append and shift mantissa
    d = d(dnonzero);
    d = min(d,INTLAB_LONG_PRECISION);
    n = size(d,1);
    precA = size(A.mantissa,2);
    precRes = max(d) + precA;
    indexRes = zeros(precRes,n);
    Amant = indexRes;
    index = precRes*(0:n-1)' + d + 1;
    indexRes(index) = 1;
    if precA~=1
      indexRes(index+precA-1) = -1;
      indexRes = logical( cumsum(indexRes,1) );
      indexRes(index+precA-1) = 1;
    else
      indexRes = logical(indexRes);
    end
    Amant(indexRes) = A.mantissa(dnonzero,:)';
    A.mantissa = [ A.mantissa zeros(size(A.mantissa,1),precRes-precA) ];
    A.mantissa(dnonzero,:) = Amant';

    % adjust exponent
    A.exponent(dnonzero) = A.exponent(dnonzero) + d;

    % adjust error and mantissa
    if precRes>INTLAB_LONG_PRECISION+1
      if INTLAB_LONG_ERROR
        A.error = errorupdate( 1 , A.error , 0 , ...
                     1 , any( A.mantissa(:,precRes+1:end),2 ) , ...
                     A.exponent-INTLAB_LONG_PRECISION-1 );
      end
      A.mantissa = A.mantissa(:,1:INTLAB_LONG_PRECISION+1);
    end

  end
