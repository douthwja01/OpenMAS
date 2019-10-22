function a = abss(a)
%ABSS         Implements  abs(a)  for intervals, result real
%
%   c = abss(a)
%
%On return, abs(alpha) <= c for all alpha in a
%Obsolete, replaced by mag (thanks to Arnold for better naming).
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 07/08/02     S.M. Rump  changed from abs->abss
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/18/08     S.M. Rump  obsolete
%

  a = mag(a);

