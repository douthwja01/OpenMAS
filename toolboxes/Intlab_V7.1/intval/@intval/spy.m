function spy(A)
%SPY          Spy sparse interval matrix
%
%   spy(A)
%

% written  10/24/99     S.M. Rump
% modified 08/07/02     S.M. Rump  abs to abss
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  spy(mag(A))
