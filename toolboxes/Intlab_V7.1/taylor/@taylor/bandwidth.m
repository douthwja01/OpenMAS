function [p,q] = bandwidth(A)
%BANDWIDTH    Upper and lower bandwidth of matrix A
%
%   A   = 0  for  i-j > p  or j-i > q
%    ij
%
%    [p,q] = bandwidth(A)
%

% written  05/21/09     S.M. Rump
%

  [p,q] = bandwidth(reshape(A.t(1,:),A.size));
