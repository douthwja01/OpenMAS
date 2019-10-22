function C = mtimes(A,B)
%MTIMES       (Entrywise) long multiplication  A * B
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 11/18/04     S.M. Rump  omit redundant code
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 06/04/02     S.M. Rump  splitting corrected (thanks to Dirk Poot)
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
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
  
  A = long(A);
  B = long(B);

  % check length of mantissa to assure components of convolution are
  % beta-digits:  pmax = floor(bitmax/beta/(beta-1));
  nA = size(A.mantissa,1);
  nB = size(B.mantissa,1);
  precA = size(A.mantissa,2);
  precB = size(B.mantissa,2);
%   pmax = floor(bitmax/INTLAB_LONG_BETA/(INTLAB_LONG_BETA-1));
pmax = floor(flintmax/INTLAB_LONG_BETA/(INTLAB_LONG_BETA-1));
  if ( precA>pmax ) & ( precB>pmax )
    % split factor
    if precA<precB
      A1 = A;
      A1.mantissa = A1.mantissa(:,1:pmax);
      if INTLAB_LONG_ERROR
        A1.error.exp = A1.error.exp -precA+pmax;
      end
      A2 = A;
      A2.mantissa = A2.mantissa(:,pmax+1:end);
      A2.exponent = A2.exponent - pmax;
      A2 = omitleadingzeros(A2);
      C = A1 * B + A2 * B;
    else
      B1 = B;
      B1.mantissa = B1.mantissa(:,1:pmax);
      if INTLAB_LONG_ERROR
        B1.error.exp = B1.error.exp -precB+pmax;
      end
      B2 = B;
      B2.mantissa = B2.mantissa(:,pmax+1:end);
      B2.exponent = B2.exponent - pmax;
      B2 = omitleadingzeros(B2);
      C = A * B1 + A * B2;
    end  
    if rndold
      setround(rndold)
    end
    return
  end

  % Compute sign and exponent
  C.sign = A.sign .* B.sign;
  C.exponent = A.exponent + B.exponent - 1;

  % Compute mantissa
  if nA==1
    C.mantissa = conv2(B.mantissa,A.mantissa);
  elseif nB==1
    C.mantissa = conv2(A.mantissa,B.mantissa);
  else
    if nA~=nB
      error('incompatible operands for long multiplication')
    end
    C.mantissa = zeros(nA,precA+precB);
    if precA>=precB
      Ones = ones(1,precA);
      for i=1:precB
        C.mantissa(:,i:i+precA-1) = C.mantissa(:,i:i+precA-1) + ...
           A.mantissa .* ( B.mantissa(:,i)*Ones );
      end
    else
      Ones = ones(1,precB);
      for i=1:precA
        C.mantissa(:,i:i+precB-1) = C.mantissa(:,i:i+precB-1) + ...
           B.mantissa .* ( A.mantissa(:,i)*Ones );
      end
    end
  end

  % Compute error
  if INTLAB_LONG_ERROR
    Amant = A.mantissa(:,1);
    if precA>1
      Amant = Amant + ( A.mantissa(:,2)+1 )/   INTLAB_LONG_BETA;
    end
    Bmant = B.mantissa(:,1);
    if precB>1
      Bmant = Bmant + ( B.mantissa(:,2)+1 )/   INTLAB_LONG_BETA;
    end
    C.error = errorupdate( Amant , B.error , A.exponent-1 , ...
                           Bmant , A.error , B.exponent-1 , ...
                           A.error.mant , B.error , A.error.exp );
  end

  % Normalize mantissa
  C = normalize(C);
  C = normalizefirst(C);

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
