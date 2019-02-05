function r = randomc(varargin)
%RANDOMC      Complex random numbers in M+i*M with M=[min,max], default [-1,+1]
%
% Calling conventions as function  random
%

% written  03/18/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if length(varargin)==0
    r = random(NaN) + sqrt(-1)*random(NaN);
  else
    r = random(NaN, varargin) + sqrt(-1)*random(NaN,varargin);
  end
  
  if rndold
    setround(rndold)
  end
