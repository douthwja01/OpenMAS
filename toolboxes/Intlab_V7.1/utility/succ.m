function r = succ(a,k)
%SUCC         successor of  a  in floating point, real or complex
%
%   r = succ(a)
%
%On return, r consists of the smallest floating point numbers with a < r
%  Elementwise for vectors and matrices,
%  with respect to partial ordering for complex entries
%
%Similarly,
%
%   r = succ(a,k)
%
%produces the  k-th successor. Negative k yields k-th predecessor.
%

% written  10/16/98     S.M. Rump
% modified 11/16/98     S.M. Rump  second parameter added, eta
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 07/08/02     S.M. Rump  allow negative k
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/03/05     S.M. Rump  comment corrected
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

  INTLAB_INTVAL_ETA = realmin*eps;

  if nargin==1
    k = 1;
  end

  r = a;
  if k<0
    setround(-1)
    for i=1:(-k)
      r = r - INTLAB_INTVAL_ETA;
    end
  else
    setround(1)
    for i=1:k
      r = r + INTLAB_INTVAL_ETA;
    end
  end
    
  setround(rndold)
