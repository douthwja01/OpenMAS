function C = rad(A)
%RAD          Radius of long (only with error term)
%
%  C = rad(A)
%

% written  11/06/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 09/20/04     S.M. Rump  output number rather than interval  
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

  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  if ~INTLAB_LONG_ERROR
    error('rad called for long without error term')
  end

  C = sup( ( sup(A) - inf(A) ) / 2 ); 
  
  if rndold
    setround(rndold)
  end
