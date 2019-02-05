function n = dim(A)
%DIM          Dimension of a square matrix
%
%    n = dim(A)
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ( length(A.size)>2 ) | ( A.size(1)~=A.size(2) )
    error('function dim called with non-square matrix')
  end;

  n = A.size(1);
