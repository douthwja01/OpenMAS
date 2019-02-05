function a = trace(a)
%TRACE        Implements  trace(a)  for slopes
%
%   c = trace(a)
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = sum( diag( a ) );
