function C = normalizefirst(C)
%NORMALIZEFIRST  Normalize first digit of long number
%
%  C = normalizefirst(C)
%
%For C.mantissa being nonnegative integers with  0 <= C.mantissa < 2^53
%with all but possibly first digit normalized to
%  0 <= C.mantissa < beta    except first digit,
%leading digit is normalized and possibly length increased
%leading digit is nonzero unless entire number is zero
%exponent and error are adapted
%
%Input C need not be of type long but only has to consist of struct fields
%

% written  12/30/98     S.M. Rump
% modified 11/19/04     S.M. Rump  error update corrected
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');

  index = ( C.mantissa(:,1)>=INTLAB_LONG_BETA );
  if any( index )
    if size(C.mantissa,2)<INTLAB_LONG_PRECISION
      C.mantissa = [ C.mantissa zeros(size(C.mantissa,1),1) ];
      C.mantissa(index,:) = [ zeros(sum(index),1) C.mantissa(index,1:end-1) ];
      C.exponent(index) = C.exponent(index) + 1;
    else
      if INTLAB_LONG_ERROR
        [ C.error.mant(index) , C.error.exp(index) ] = ...
             errorupdate( 1 , C.error.mant(index) , C.error.exp(index) , ...
                         1 , C.mantissa(index,end) , ...
                         C.exponent(index)-INTLAB_LONG_PRECISION );
      end
      C.mantissa(index,:) = [ zeros(sum(index),1) C.mantissa(index,1:end-1) ];
      C.exponent(index) = C.exponent(index) + 1;
    end
    C.mantissa(index,1) = floor( C.mantissa(index,2)/INTLAB_LONG_BETA );
    C.mantissa(index,2) = C.mantissa(index,2) - ...
                          C.mantissa(index,1)*INTLAB_LONG_BETA;
  end
