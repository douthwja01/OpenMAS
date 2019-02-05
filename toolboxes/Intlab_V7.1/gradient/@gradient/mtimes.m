function r = mtimes(a,b)
%MTIMES       Gradient multiplication  a * b
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    case distinction for gradient / non-gradient input
%                                    complete redesign
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/02/05     S.M. Rump  improved performance (thanks to Joerg Kubitz, Hannover)
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  if ( length(size(a))>2 ) | ( length(size(b))>2 )
    error('gradient multiplication * only for 2-dimensional arrays')
  end
  [m k] = size(a);              % first dimension
  [k1 n] = size(b);             % second dimension
  if m*k==1                     % one factor scalar
    r = a .* b;
    return
  elseif k1*n==1                % one factor scalar
    r = b .* a;
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  % both factors arrays
  if k~=k1
    error('inner dimensions do not match')
  end
  N = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  if ~isa(a,'gradient')      % non-gradient array * gradient array

    INTLAB_GRADIENT_SPARSE = getappdata(0,'INTLAB_GRADIENT_SPARSE');
    if N>=INTLAB_GRADIENT_SPARSE
      a = sparse(a);
    end
    r.x = a * b.x;
    e = b.dx.';
    e = reshape(reshape(e(:,reshape(1:k*n,k,n)'),N*n,k)*(a.'),N,m*n);
    r.dx = e(:,reshape(1:n*m,n,m)').';

  elseif ~isa(b,'gradient')          % gradient array * non-gradient array

    INTLAB_GRADIENT_SPARSE = getappdata(0,'INTLAB_GRADIENT_SPARSE');
    if N>=INTLAB_GRADIENT_SPARSE
      b = sparse(b);
    end
    r.x = a.x * b;
    r.dx = reshape(reshape(a.dx.',N*m,k)*b,N,m*n).';

  else                                  % gradient array * gradient array

    r.x = a.x * b.x;
    e = b.dx.';
    dims = m*n;
    e = reshape(reshape(e(:,reshape(1:k*n,k,n)'),N*n,k)*(a.x).',N,dims);
    r.dx = ( e(:,reshape(1:n*m,n,m)') + reshape(reshape(a.dx.',N*m,k)*(b.x),N,dims) ).';

  end

  if m*n==1                     % avoid Matlab sparse 1x1 bug
    r.x = full(r.x);
    if N==1
      r.dx = full(r.dx);
    end
  end
  r = class(r,'gradient');

  if rndold
    setround(rndold)
  end
