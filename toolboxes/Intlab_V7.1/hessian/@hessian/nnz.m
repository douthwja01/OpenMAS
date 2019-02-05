function nz = nnz(a)
%NNZ          Number of nonzero elements in hessian
%
%  nz = nnz(a)
%
%Result is (1,3) vector [ nnz(a.x) nnz(a.dx) nnz(a.hx) ].
%

% written  04/04/04     S.M. Rump 
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  index = reshape(1:N^2,N,N)';
  ahx = a.hx + a.hx(index(:),:);     % Hessian is .hx + transpose(.hx)
  nz = [ nnz(a.x) nnz(a.dx) nnz(ahx) ];
  