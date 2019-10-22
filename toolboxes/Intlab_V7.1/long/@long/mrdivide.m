function C = mrdivide(A,B)
%MRDIVIDE     Long division A / B
%

% written  12/30/98     S.M. Rump
% modified 10/22/99     S.M. Rump  improvement of error if 1/B exact in few digits
% modfied  02/09/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 01/29/10     S.M. Rump  midpoint corrected, thanks to Nozomu Matsuda and Nobito Yamamoto
% modified 06/09/10     S.M. Rump  division w/o error term, thanks to Nozomu Matsuda
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
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
  
  A = long(A);

  if isa(B,'double') & all( B==round(B) ) & all( abs(B)<INTLAB_LONG_BETA )

    % check zero denominator
    if any( B==0 )
      error('long division by zero')
    end

    % Division by integer
    C = A;
    C.sign = C.sign .* sign(B);
    B = abs(B);

    % Compute mantissa
    precA = size(A.mantissa,2);
    C.mantissa(:,1) = floor( A.mantissa(:,1)./B );
    for i=2:precA
      A.mantissa(:,i) = A.mantissa(:,i) + ...
         ( A.mantissa(:,i-1) - B.*C.mantissa(:,i-1) )*INTLAB_LONG_BETA ;
      C.mantissa(:,i) = floor( A.mantissa(:,i)./B );
    end

    % If remainder nonzero, continue division
    rem = A.mantissa(:,precA) - B.*C.mantissa(:,precA);
    if any( rem )
      C.mantissa = ...
         [ C.mantissa zeros(size(C.mantissa,1),INTLAB_LONG_PRECISION-precA+1) ];
      i = precA;
      while any( rem ) & ( i<INTLAB_LONG_PRECISION+1 )
        i = i+1;
        num = rem*INTLAB_LONG_BETA;
        C.mantissa(:,i) = floor( num./B );
        rem = num - B.*C.mantissa(:,i);
      end
    end

    % Omit trailing zeros
    if i<INTLAB_LONG_PRECISION
      C.mantissa = C.mantissa(:,1:i);
    end

    % Compute error
    if INTLAB_LONG_ERROR
      precC = size(C.mantissa,2);
      C.error = errorupdate( -B , C.error , 0 , 1 , rem~=0 , C.exponent-precC );
    end

    % Omit leading zeros
    C = omitleadingzeros(C);

    % Normalize result
    C = normalize(C);
    if INTLAB_LONG_ERROR
      C.error = errornormalize(C.error);
    end

  else   % long division

    % denominator with one extra digit precision
    INTLAB_LONG_PRECISION = INTLAB_LONG_PRECISION + 1;
    B = long(B);
    if INTLAB_LONG_ERROR
      Bmid = mid(B);
      % denominator offset rad(B)^2/mid(B) in double
      N = long2intval(rad(B)).^2./long2intval(Bmid);
      C = addlongerror(Bmid - mid(N),rad(N));     % to be inverted
    else
      Bmid = B;
      C = Bmid;
    end
    precC = size(C.mantissa,2);

    % check zero denominators
    if any( C.mantissa(:,1)==0 )
      error('long division by zero')
    end

    % approximation of inverse for Newton iteration, approx. 52 bits accuracy
    Cinv = C;
    Cinv.exponent = 0;
    Cinv = long(1./long2dble(Cinv));    % inverse of mid(C)
    Cinv.exponent = Cinv.exponent - C.exponent;

    % approximate inverse without error term
    err = longinit('errorterm',0);
    longinit('WithoutErrorTerm',0);

    while 1
      Res = 2 - Cinv*C;
      % check Res==1
      if all( Res.mantissa(:,1)==1 ) & all(all( Res.mantissa(:,2:end)==0 ))
        break
      end
      Cinv = Cinv * Res;
    end

    % reset precision and error term, Cinv~1/C
    longinit(err,0);
    INTLAB_LONG_PRECISION = INTLAB_LONG_PRECISION - 1;

    if INTLAB_LONG_ERROR
      % calculate error
      % error by inversion:  Res_/C  with
      %   Res_ < beta^(1-precRes),   1/C >= 1/C_1*beta^(1-C.exponent)
      % even Res < beta^(-INTLAB_LONG_PRECISION) because Newton iteration
      %   stopped with Res==1 in precision INTLAB_LONG_PRECISION+1; this
      %   is important if 1/C is exactly representable in fewer than
      %   INTLAB_LONG_PRECISION+1 beta-digits (thanks to Kurt Zehetleitner)
      precCinv = size(Cinv.mantissa,2);
      Cmant = C.mantissa(:,1);
      if precC>1
        Cmant = Cmant + C.mantissa(:,2)/INTLAB_LONG_BETA;
      end
      Cinv.error = ...
           errorupdate( -Cmant , 1 , 1-C.exponent-INTLAB_LONG_PRECISION );

      % error introduced by B:  B.error/(sqr(B)-sqr(B.error))
      Bmant = B.mantissa(:,1);
      if size(B.mantissa,2)>1
        Bmant = Bmant + B.mantissa(:,2)/INTLAB_LONG_BETA;
      end
      setround(1)
      N = ( B.error.mant .* INTLAB_LONG_BETA.^(B.error.exp-B.exponent) ) .^ 2;
      setround(-1)
      N = ( Bmant/INTLAB_LONG_BETA ) .^ 2 - N;
      % true denominator  N*beta^(2*B.exponent)

      % check zero denominator (including error)
      if any( N<=0 )
        error('long division by zero')
      end

      % first part:   B.error/(sqr(B)-sqr(B.error))
      Cinv.error = errorupdate( 1 , Cinv.error , 0 , ...
                               -N , B.error , -2*B.exponent );
      Cinv.error = errornormalize(Cinv.error);
    end

    C = A*Cinv;

  end
  
  setround(rndold)
