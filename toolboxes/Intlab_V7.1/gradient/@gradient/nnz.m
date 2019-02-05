function nz = nnz(a)
%NNZ          Number of nonzero elements in gradient
%
%  nz = nnz(a)
%
%Result is (1,3) vector [ nnz(a.x) nnz(a.dx) ].
%

% written  04/08/04     S.M. Rump 
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  nz = [ nnz(a.x) nnz(a.dx) ];
  