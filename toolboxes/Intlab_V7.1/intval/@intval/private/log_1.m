function x = log_1(x,rnd)
%LOG_         Rigorous calculation of  log(1+x)  for 0<=x<=1 according to rnd
%
%   y = log_1(x)
%
%Internal function
%

% written  12/30/98     S.M. Rump
% modified 08/31/98     S.M. Rump  improved accuracy
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_LOG = getappdata(0,'INTLAB_STDFCTS_LOG');

  setround(0)
  xs = pow2( floor(x*2^13) , -13 );    % first 13 bits
  log1xs = log(1+xs);                  % 1+xs exactly representable in 14 bits

  setround(rnd)
  d = ( x - xs ) ./ (1+xs);            % 0 <= d < 2^-13

  % log(1+x) = log( (1+xs) * (1+d) ) ,   0 <= err <= d^5/5 < 4.45e-17*d
  if rnd==-1
    log1d = ((( (-d)/4 + 1/3 ).*d - 0.5 ).*d).*d + d;
    x = log1xs + ( log1d + (-INTLAB_STDFCTS_LOG.EPS)*abs(log1xs) );
  else
    log1d = (((( d/5 - .25 ).*d + 1/3 ).*d - 0.5 ).*d).*d + d;
    x = log1xs + ( log1d + INTLAB_STDFCTS_LOG.EPS*abs(log1xs) );
  end

  setround(0)
