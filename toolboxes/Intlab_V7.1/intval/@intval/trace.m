function c = trace(a)
%TRACE        Implements  trace(a)  for interval matrices
%
%   c = trace(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  c = sum( diag( a ) );
