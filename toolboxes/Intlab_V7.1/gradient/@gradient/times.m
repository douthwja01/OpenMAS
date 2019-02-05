function r = times(a,b)
%TIMES        Gradient multiplication  a .* b
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  array/scalar operations
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
%                                    complete redesign
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/05/05     S.M. Rump  improved performance (thanks to J. Kubitz)
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

  N = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  if ~isa(a,'gradient')
    m = prod(size(a));
    if m==1                     % non-gradient scalar .* gradient
      r.x = a * b.x;
      r.dx = a * b.dx;
    else                        % non-gradient array .* gradient
      n = prod(size(b.x));
      if n==1                   % non-gradient array .* gradient scalar
        r.x = a * b.x;
        ax = a(:);
        if issparse(b.dx)
          ax = sparse(ax);
        end
        r.dx = ax * b.dx;
      else                      % non-gradient array .* gradient array
        if ~isequal(size(a),size(b))
          error('gradient .* : dimensions not compatible')
        end
        r.x = a .* b.x;
        if issparse(b.dx)
          [ib,jb,sb] = find(b.dx);
          r.dx = sparse(ib,jb,reshape(a(ib),size(sb)).*sb,n,N);
        else
          r.dx = b.dx .* repmat(a(:),1,N);
        end          
      end
    end
  elseif ~isa(b,'gradient')      % gradient times non-gradient
    r = b .* a;
    if rndold
      setround(rndold)
    end
    return
  else                          % both factors gradient
    m = prod(size(a.x));
    n = prod(size(b.x));
    sparse_ = issparse(a.dx) | issparse(b.dx);
    if m==1                     % scalar gradient .* gradient
      if n==1                   % scalar gradient .* scalar gradient
        r.x = a.x * b.x;
        r.dx = a.dx * b.x + a.x * b.dx;
      else                      % scalar gradient .* array gradient
        r.x = a.x * b.x;
        if sparse_
          r.dx = sparse(b.x(:)) * a.dx + b.dx * a.x;
        else
          r.dx = b.x(:) * a.dx + b.dx * a.x;
        end
      end
    else                        % array gradient .* gradient
      if n==1                   % array gradient .* scalar gradient
        r.x = a.x * b.x;
        if sparse_
          r.dx = sparse(a.x(:)) * b.dx + a.dx * b.x;
        else
          r.dx = a.x(:) * b.dx + a.dx * b.x;
        end
      else                      % array gradient .* array gradient
        if ~isequal(size(a),size(b))
          error('dimensions not compatible for gradient .*');
        end
        r.x = a.x .* b.x;
        if issparse(a.dx) | issparse(b.dx)
          r.dx = sparse(1:n,1:n,b.x(:))*a.dx + sparse(1:n,1:n,a.x(:))*b.dx;
        else
          r.dx = repmat(b.x(:),1,N) .* a.dx + b.dx .* repmat(a.x(:),1,N);
        end
      end
    end
  end

  r = class(r,'gradient');
  
  if rndold
    setround(rndold)
  end
