function a = sparse(a)
%SPARSE       Convert hessian first and second derivative to sparse
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/06/05     S.M. Rump  .x part also sparse
% modified 10/03/12     S.M. Rump  sparse up to 2 dimensions
%

  if length(size(a))>2
    error('sparse arrays only up to 2 dimensions')
  end
  
  % avoid Matlab 6.5f bug:
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if prod(size(a.x))~=1
    a.x = sparse(a.x);
  end
  if prod(size(a.dx))~=1
    a.dx = sparse(a.dx);
  end
  if prod(size(a.hx))~=1
    a.hx = sparse(a.hx);
  end
