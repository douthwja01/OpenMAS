function m = shiftright(m,r)
%SHIFTRIGHT   Shift array m right by r bits
%
%   m = shiftright(m,r)
%
%last r bits of m(:,end) are zero, otherwise lost
%

% written  12/30/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');

  % array m divided by 2^r
  q = m .* ( 2.^(-r)*ones(1,size(m,2)) );

  % array m shifted right by r bits (last r bits lost)
  floorq = floor(q);

  % first component of m
  m(:,1) = floorq(:,1);

  % shifted array plus last r bits of previous array component
  m(:,2:end) = ( q(:,1:end-1) - floorq(:,1:end-1) ) * INTLAB_LONG_BETA + ...
               floorq(:,2:end);
