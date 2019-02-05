function res = CompMat(a);
%COMPMAT      Ostrowski's comparison matrix
%
%   res = CompMat(a)
%

% written   8/12/94     S.M. Rump
% modified 09/22/02     S.M. Rump  check square matrix and interval matrices
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  [m n] = size(a);
  if m~=n
    error('Comparison matrix only for square matrices')
  end
  
  res = -mag(a);
  res(1:n+1:n*n) = mig(diag(a));
