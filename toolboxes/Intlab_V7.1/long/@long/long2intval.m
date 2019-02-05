function c = long2intval(C)
%LONG2INTVAL  Conversion long to intval (with correct rounding)
%
%  c = long2intval(C)
%

% written  12/30/98     S.M. Rump
% modfied  02/09/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 07/31/05     S.M. Rump  header corrected
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
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

  % extra treatment of zero component
  indexzero = ( C.exponent==-inf );

  n = size(C.mantissa,1);
  % take maximum first 80 bits
  precC = size(C.mantissa,2);
  p = min( precC , ceil(80/INTLAB_LONG_LOGBETA)+1 );
  if INTLAB_LONG_ERROR
    if p<precC
      C.error = errorupdate( 1 , C.error , 0 , ...
                  1 , any(C.mantissa(:,p+1:precC)~=0,2) , C.exponent-precC );
    end
    Cerr = C.error.mant .* INTLAB_LONG_BETA.^(C.error.exp-C.exponent+p);
  else
    Cerr = 0;
  end

  factor = ( INTLAB_LONG_BETA.^(-p:-1) )';
  exppos = ( C.exponent>=0 );
  expneg = ~exppos;
  E = floor( 600/INTLAB_LONG_LOGBETA );
  F = INTLAB_LONG_BETA ^ E;

  Cmantp = C.mantissa(:,p);
  setround(-1)
  C.mantissa(:,p) = Cmantp - Cerr;
  cinf = C.mantissa(:,p:-1:1) * factor ;
  if any(exppos)
    cinf(exppos) = ( cinf(exppos)*F ) .* ...
                      ( INTLAB_LONG_BETA.^(C.exponent(exppos)-E) );
  end
  if any(expneg)
    cinf(expneg) = ( cinf(expneg)/F ) .* ...
                      ( INTLAB_LONG_BETA.^(C.exponent(expneg)+E) );
  end

  setround(1)
  C.mantissa(:,p) = Cmantp + Cerr;
  csup = C.mantissa(:,p:-1:1) * factor ;
  if any(exppos)
    csup(exppos) = ( csup(exppos)*F ) .* ...
                      ( INTLAB_LONG_BETA.^(C.exponent(exppos)-E) );
  end
  if any(expneg)
    csup(expneg) = ( csup(expneg)/F ) .* ...
                      ( INTLAB_LONG_BETA.^(C.exponent(expneg)+E) );
  end
  if INTLAB_LONG_ERROR & any(indexzero)
    csup(indexzero) = C.error.mant(indexzero) .* ...
                        INTLAB_LONG_BETA.^C.error.exp(indexzero);
    cinf(indexzero) = - csup(indexzero);
  end

  c = hull( C.sign.*cinf , C.sign.*csup );
  
  setround(rndold)
