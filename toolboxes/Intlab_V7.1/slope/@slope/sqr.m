function u = sqr(a)
%SQR          Slope square  sqr(a)
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  u = a;

  u.r = sqr(a.r);
  indexc = 1:INTLAB_SLOPE.NUMVAR;
  indexr = 2:INTLAB_SLOPE.NUMVAR+1;
  u.s = a.s .* (a.r(:,indexc)+a.r(:,indexr));
  
  if rndold
    setround(rndold)
  end
