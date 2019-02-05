function a = sparse(a)
%SPARSE       Convert Taylor to sparse
%

% written  05/21/09     S.M. Rump
%

  a.t = sparse(a.t);
  