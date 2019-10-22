function res = subrealmin
%SUBREALMIN   Smallest unnormalized positive floating-point number
%
%   res = subrealmin
%

% written  08/07/10     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  res = realmin*eps;
