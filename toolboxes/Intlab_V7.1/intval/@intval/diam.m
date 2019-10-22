function c = diam(a)
%DIAM         implements  diam(a)
%
%   c = diam(a)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  treats infinity intervals
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 04/06/05     S.M. Rump  fast check for rounding to nearest
%

  if a.complex
    c = 2*a.rad;
  else
    e = 1e-30;
    if 1+e==1-e                   % fast check for rounding to nearest
      rndold = 0;
    else
      rndold = getround;
    end
    setround(1)
    c = a.sup - a.inf;
    setround(rndold)              % set rounding to previous value
  end
