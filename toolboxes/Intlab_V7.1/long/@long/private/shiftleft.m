function m = shiftleft(m,r)
%SHIFTLEFT    Shift array m left by r bits
%
%   m = shiftleft(m,r)
%
%first r bits of m(:,1) are zero, otherwise lost
%

% written  12/30/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');

  Ones = ones(1,size(m,2));

  % first r bits of m
  factor = 2.^(r-INTLAB_LONG_LOGBETA)*Ones ;
  ms = floor( m .* factor );

  % array m shifted left by r bits, last r bits zero
  m =  ( m - ms ./ factor ) .* ( 2.^r*Ones );

  % add last r bits
  m(:,1:end-1) = m(:,1:end-1) + ms(:,2:end);
