function n = dim(A)
%DIM          Dimension of a square matrix
%
%    n = dim(A)
%

% written  05/21/09     S.M. Rump
%

  n = A.size(1);
  if ( length(A.size)~=2 ) | ( n ~= A.size(2) )
    error('function dim called with non-square matrix')
  end;
