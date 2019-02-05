function C = addlongerror(C,err,exponent)
%ADDLONGERROR Add error to long number
%
%  C = addlongerror(A,err,exponent)
%
%Error term of C is updated by err. Third argument optional, default 0.
%If exponent is specified, error term is updated by err*10^exponent
%
%Input C, err and exponent may be column vectors
%

% written  12/30/98     S.M. Rump
% modfied  02/09/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 01/29/10     S.M. Rump  only warning if no interval arithmetic 
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  
  if ~INTLAB_LONG_ERROR
    warning('long arithmetic changed to interval arithmetic, see longinit')
    INTLAB_LONG_ERROR = 1;
  end

  if any( err<0 )
    error('negative error specified in addlongerror')
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==3
    setround(1)
    % 10^exponent = beta^(Ernd+Efrac)
    E = exponent / log10(INTLAB_LONG_BETA);
    Ernd = round(E);
    Efrac = E - Ernd;
    err = err .* INTLAB_LONG_BETA.^Efrac;
    setround(0)                         % set rounding to nearest
  else
    Ernd = 0;
  end
  C.error = errorupdate( 1 , C.error , 0  , 1 , err , Ernd );
  C.error = errornormalize(C.error);

  C = normalize(C);
  
  setround(rndold)
