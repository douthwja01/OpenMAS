function ur = rangeimprove(u)
%RANGEIMPROVE Possible sharpening of range using slope and intersection
%
%Internal function

%Input    structure u with u.r and u.s
%Output   improved range ur by intersection

% written  12/06/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  M = [ u.r(:,1) u.s .* repmat(INTLAB_SLOPE.Xxs',size(u.s,1),1) ];
  setround(-1)
  urinf = cumsum(inf(M),2);
  setround(1)
  ursup = cumsum(sup(M),2);
  ur = intersect( u.r , infsup(urinf,ursup) );
