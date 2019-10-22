function r = issparse(c)
%ISSPARSE     Returns 1 if c is sparse
%
%  r = issparse(c)
%

% written  05/21/09     S.M. Rump
%

  r = issparse(c.t);
