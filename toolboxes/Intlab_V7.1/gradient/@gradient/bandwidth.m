function [p,q] = bandwidth(A)
%BANDWIDTH    Upper and lower bandwidth of matrix A
%
%   A   = 0  for  i-j > p  or j-i > q
%    ij
%
%    [p,q] = bandwidth(A)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  [p,q] = bandwidth(A.x);
