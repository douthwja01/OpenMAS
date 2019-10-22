function c = cintval(a,r)
%CINTVAL      Initialization of complex interval
%
%   c = cintval(a,r)
%
%Result c computed such that <a,r> is enclosed in interval  c.
%  Input argument r is optional, default 0.
%
%Output interval is complex, i.e.  c = { z in C :  |z-a| <= r }
%For complex a, this is the same as midrad(a,r).
%

% written  11/23/98     S.M. Rump
% modified 06/24/99     S.M. Rump  check sparsity, multi-dimensional arrays
% modified 09/02/00     S.M. Rump  replaces midradcmplx
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/18/08     S.M. Rump  radius zero for sparse input
%

  if nargin==1
    r = 0;
  end

  index = ( r<0 );
  if ~isreal(r) | any(index(:))
    error('invalid radius in call of cintval')
  end

  if ( ~isequal(r,0) ) & ( issparse(a)~=issparse(r) )
    error('midpoint and radius must both be or none of them be sparse')
  end

  c = intval(a,r,'midrad');
