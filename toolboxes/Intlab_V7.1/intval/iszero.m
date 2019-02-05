function r = iszero(a)
%ISZERO       Returns array of 1's for zero components (for completeness)
%
%   r = iszero(a)
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/12/05     S.M. Rump  typo corrected
%

  r = ( a==0 );
