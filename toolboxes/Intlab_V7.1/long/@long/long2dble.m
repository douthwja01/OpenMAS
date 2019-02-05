function c = long2dble(C)
%LONG2DBLE    Conversion long to double (rounding to nearest), w/o eror term
%
%  c = long2dble(C)
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
  
  n = size(C.mantissa,1);
  % take maximum first 80 bits
  p = min( size(C.mantissa,2) , ceil(80/INTLAB_LONG_LOGBETA)+1 );

  factor = ( INTLAB_LONG_BETA.^(-p:-1) )';
  exppos = ( C.exponent>=0 );
  expneg = ~exppos;
  E = floor( 600/INTLAB_LONG_LOGBETA );
  F = INTLAB_LONG_BETA ^ E;

  % calculate mantissa for small and large exponents
  c = C.mantissa(:,p:-1:1) * factor ;
  if any(exppos)
    c(exppos) = ( c(exppos)*F ) .* ( INTLAB_LONG_BETA.^(C.exponent(exppos)-E) );
  end
  if any(expneg)
    c(expneg) = ( c(expneg)/F ) .* ( INTLAB_LONG_BETA.^(C.exponent(expneg)+E) );
  end

  % get the sign
  c = C.sign .* c;
  
  if rndold
    setround(rndold)
  end
