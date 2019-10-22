function c = rad(a)
%RAD          Implements  rad(a)  for intervals
%
%   c = rad(a)
%
% mid(a) and rad(a) computed such that
%    alpha  in  < mid(a) , rad(a) >  for all alpha in a      (*)
%
%As has been pointed out by G. Mayer, Rostock, this implies a peculiarity.
%For an interval X with bound differing by one bit, the result of rad(X)
%and diam(X) is the same. For example, 
%
%   x = infsup(1,1+eps), rx = rad(x), dx = diam(x)
%
%yields
% 
% intval x = 
%     1.0000
% rx =
%   2.2204e-016
% dx =
%   2.2204e-016
%
%to assure (*). So for comparing the quality of interval inclusions I recommend
%  to use diam.
%

% written  10/16/98     S.M. Rump
% modified 06/22/99     S.M. Rump  for sparse matrices
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 06/29/05     S.M. Rump  comment rad/diam added
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 09/10/07     S.M. Rump  performance, huge arrays
% modified 07/23/09     S.M. Rump  changed formula (thanks to Gerhard Heindl for pointing to this)
%

  if a.complex
    c = a.rad;
  else
    e = 1e-30;
    if 1+e==1-e                         % fast check for rounding to nearest
      rndold = 0;
    else
      rndold = getround;
    end
    m = mid(a);
    setround(1)
    c = max( m-a.inf , a.sup-m );
    setround(rndold)
  end
