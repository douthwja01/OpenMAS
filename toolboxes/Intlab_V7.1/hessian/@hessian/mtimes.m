function r = mtimes(a,b)
%MTIMES       Hessian multiplication  a * b
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/05/05     S.M. Rump  improved performance
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/10/08     S.M. Rump  make sure r.hx remainse sparse
% modified 08/26/12     S.M. Rump  global variables removed
%

  if ( length(size(a))>2 ) | ( length(size(b))>2 )
    error('hessian multiplication * only for 2-dimensional arrays')
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
  
  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');
  
  if ~isa(a,'hessian')          % non-hessian array * hessian array
    
    INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
    if N>=INTLAB_HESSIAN_SPARSE
      a = sparse(a);
    end
    r.x = a * b.x;
    dims = m*n;
    r.dx = reshape(reshape(b.dx(:,reshape(1:k*n,k,n)'),N*n,k)*(a.'),N,dims);
    r.dx = r.dx(:,reshape(1:n*m,n,m)');
    r.hx = reshape(reshape(b.hx(:,reshape(1:k*n,k,n)'),N^2*n,k)*(a.'),N^2,dims);
    r.hx = r.hx(:,reshape(1:n*m,n,m)');

  elseif ~isa(b,'hessian')      % hessian array * non-hessian array
    
    INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
    if N>=INTLAB_HESSIAN_SPARSE
      b = sparse(b);
    end
    r.x = a.x * b;
    dims = m*n;
    r.dx = reshape(reshape(a.dx,N*m,k)*b,N,dims);
    r.hx = reshape(reshape(a.hx,N^2*m,k)*b,N^2,dims);

  else                                  % hessian array * hessian array
    
    r.x = a.x * b.x;                    % hessian array * hessian array (sparse)
    dims = m*n;
    r.dx = reshape(reshape(b.dx(:,reshape(1:k*n,k,n)'),N*n,k)*(sparse(a.x)).',N,dims);
    r.dx = r.dx(:,reshape(1:n*m,n,m)') + reshape(reshape(a.dx,N*m,k)*(sparse(b.x)),N,dims);
    e = reshape(b.dx.',k,n*N);
    r.hx = reshape(a.dx,m*N,k) * reshape(e(:,reshape(1:n*N,n,N)'),k,n*N);
    e = reshape(reshape(b.hx(:,reshape(1:k*n,k,n)'),N^2*n,k)*(sparse(a.x)).',N^2,dims);
    r.hx = reshape(reshape(r.hx(:,reshape(1:n*N,N,n)'),m*n*N,N).',N^2,m*n) + ...
      e(:,reshape(1:n*m,n,m)') + reshape(reshape(a.hx,N^2*m,k)*sparse(b.x),N^2,dims);

  end
  
  if m*n==1                     % avoid Matlab sparse 1x1 bug
    r.x = full(r.x);
    if N==1
      r.dx = full(r.dx);
      r.hx = full(r.hx);
    end
  end
  r = class(r,'hessian');
  
  if rndold
    setround(rndold)
  end
