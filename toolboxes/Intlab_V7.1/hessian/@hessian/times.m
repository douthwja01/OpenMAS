function r = times(a,b)
%TIMES        Hessian elementwise multiplication  a .* b
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  if ~isa(a,'hessian')
    m = prod(size(a));
    if m==1                     % non-hessian scalar .* hessian
      r.x = a * b.x;
      r.dx = a * b.dx;
      r.hx = a * b.hx;
    else                        % non-hessian array .* hessian
      n = prod(size(b.x));
      if n==1                   % non-hessian array .* hessian scalar
        r.x = a * b.x;
        ax = a(:).';
        if issparse(b.hx)
          ax = sparse(ax);
        end
        r.dx = b.dx * ax;
        r.hx = b.hx * ax;
      else                      % non-hessian array .* hessian array
        if ~isequal(size(a),size(b))
          error('hessian .* : dimensions not compatible')
        end
        r.x = a .* b.x;
        if issparse(b.hx)
          [ib,jb,sb] = find(b.dx);
          r.dx = sparse(ib,jb,reshape(a(jb),size(sb)).*sb,N,n);
          [ib,jb,sb] = find(b.hx);
          r.hx = sparse(ib,jb,reshape(a(jb),size(sb)).*sb,N^2,n);
        else
          ax = repmat(a(:).',N^2,1);
          r.dx = ax(1:N,:) .* b.dx;
          r.hx = ax .* b.hx;
        end          
      end
    end
  elseif ~isa(b,'hessian')      % hessian times non-hessian
    r = b .* a;  
    if rndold
      setround(rndold)
    end
    return
  else                          % both factors hessian
    m = prod(size(a.x));
    n = prod(size(b.x));
    sparse_ = issparse(a.hx) | issparse(b.hx);
    if m==1                     % scalar hessian .* hessian
      if n==1                   % scalar hessian .* scalar hessian
        r.x = a.x * b.x;
        r.dx = a.dx * b.x + a.x * b.dx;
        if sparse_
          r.hx = a.hx * b.x + reshape(sparse(a.dx) * sparse(b.dx.'),N^2,1) + a.x * b.hx;
        else
          r.hx = a.hx * b.x + reshape(a.dx * b.dx.',N^2,1) + a.x * b.hx;
        end
      else                      % scalar hessian .* array hessian
        r.x = a.x * b.x;
        if sparse_
          bx = sparse(b.x(:).');
          r.dx = a.dx * bx + b.dx * a.x;
          r.hx = a.hx * bx + reshape(sparse(a.dx)*sparse(b.dx(:).'),N^2,n) + a.x * b.hx ;
        else
          r.dx = a.dx * b.x(:).' + b.dx * a.x;
          r.hx = a.hx * b.x(:).' + reshape(a.dx*(b.dx(:).'),N^2,n) + a.x * b.hx ;
        end
      end
    else                        % array hessian .* hessian
      if n==1                   % array hessian .* scalar hessian
        r.x = a.x * b.x;
        if sparse_
          ax = sparse(a.x(:).');
          r.dx = b.dx * ax + a.dx * b.x;
          r.hx = b.hx * ax + reshape(sparse(b.dx)*sparse(a.dx(:).'),N^2,m) + b.x * a.hx;
        else
          r.dx = b.dx * a.x(:).' + a.dx * b.x;
          r.hx = b.hx * a.x(:).' + reshape(b.dx*(a.dx(:).'),N^2,m) + b.x * a.hx;
        end
      else                      % array hessian .* array hessian
        if ~isequal(size(a),size(b))
          error('dimensions not compatible for hessian .*');
        end
        r.x = a.x .* b.x;
        if issparse(a.hx) | issparse(b.hx)
          [ia,ja,sa] = find(a.dx);
          [ib,jb,sb] = find(b.dx);
          r.dx = sparse(ia,ja,reshape(b.x(ja),size(sa)).*sa,N,m) + ...
                 sparse(ib,jb,reshape(a.x(jb),size(sb)).*sb,N,m);
          r.hx = adx2rhx(N,m,a.dx,b.dx);
          [ia,ja,sa] = find(a.hx);
          [ib,jb,sb] = find(b.hx);
          r.hx = r.hx + sparse(ia,ja,reshape(b.x(ja),size(sa)).*sa,N^2,m) + ...
                        sparse(ib,jb,reshape(a.x(jb),size(sb)).*sb,N^2,m);
        else
          ax = repmat(a.x(:).',N^2,1);
          bx = repmat(b.x(:).',N^2,1);
          r.dx = a.dx .* bx(1:N,:) + ax(1:N,:) .* b.dx;
          index = repmat( 1:N , N , 1 );
          r.hx = a.hx .* bx + ax .* b.hx + a.dx(index,:) .* b.dx(index',:);
        end
      end
    end
  end

  r = class(r,'hessian');
  
  if rndold
    setround(rndold)
  end
