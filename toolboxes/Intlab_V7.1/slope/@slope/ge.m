function res = ge(a,b)
%GE           Implements  a >= b  for slopes, compares only a.r and b.r
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isa(a,'slope')
    asize = size(a);
    bsize = b.size;
    res = ( a(:)>=b.r(:,1) );
  elseif ~isa(b,'slope')
    asize = a.size;
    bsize = size(b);
    res = ( a.r(:,1)>=b(:) );
  else
    asize = a.size;
    bsize = b.size;
    res = ( a.r(:,1)>=b.r(:,1) );
  end

  if ~isequal(asize,bsize) & ( prod(asize)~=1 ) & ( prod(bsize)~=1 )
    error('incompatible size for slope >= ')
  end

  if prod(asize)==1
    res = reshape(res,bsize);
  else
    res = reshape(res,asize);
  end
