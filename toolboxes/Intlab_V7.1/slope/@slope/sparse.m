function a = sparse(a)
%SPARSE       Convert slope to sparse
%

% written  04/06/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/03/12     S.M. Rump  sparse up to 2 dimensions
%

  if length(a.size)>2
    error('sparse arrays only up to 2 dimensions')
  end
  a.r = sparse(a.r);
  a.s = sparse(a.s);
