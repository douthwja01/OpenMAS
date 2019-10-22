function C = normalize(C)
%NORMALIZE    Normalize long number
%
%  C = normalize(C)
%
%For C.mantissa being nonnegative integers with  0 <= C.mantissa < 2^53,
%output is normalized C such that
%  0 <= C.mantissa < beta    except first digit
%leading digit is nonzero unless entire number is zero
%exponent and error are adapted
%
%Input C need not be of type long but only has to consist of struct fields
%

% written  12/30/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  while any(any( C.mantissa(:,2:end)>=INTLAB_LONG_BETA ))
    mant = C.mantissa(:,2:end);
    q = floor( mant/INTLAB_LONG_BETA );
    C.mantissa(:,2:end) = mant - q*INTLAB_LONG_BETA;
    C.mantissa(:,1:end-1) = C.mantissa(:,1:end-1) + q;
  end

  if INTLAB_LONG_ERROR
    % shorten mantissa if error too big
    precC = size(C.mantissa,2);
    q = min( ( C.error.exp - ( C.exponent-precC ) ).* ( C.error.mant~=0 ) ) - 1;
    if q>0        % error too big, skip last digits
      q = min( q , size(C.mantissa,2)-1 );
      C.error = errorupdate( 1 , C.error , 0 , ...
                1 , any( C.mantissa(:,end-q+1:end)~=0 , 2 ) , ...
                C.exponent-precC+q );
      C.mantissa = C.mantissa(:,1:end-q);
    end
  end
