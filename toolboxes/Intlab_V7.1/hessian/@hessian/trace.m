function a = trace(a)
%TRACE        Implements  trace(a)  for hessians
%
%   c = trace(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = sum( diag( a ) );
