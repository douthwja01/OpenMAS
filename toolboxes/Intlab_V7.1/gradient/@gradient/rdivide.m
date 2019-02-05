function r = rdivide(a,b)
%RDIVIDE      Gradient elementwise right division a .\ b
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  array/scalar operations
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    complete redesign
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/05/05     S.M. Rump  improved performance
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  if prod(size(b))==1               % scalar denominator
    r = a / b;
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  scalar_a = ( prod(size(a))==1 );
  if ( ~scalar_a ) & ( ~isequal(size(a),size(b)) )
    error('dimensions do not match in gradient elementwise division')
  end
  
  N = getappdata(0,'INTLAB_GRADIENT_NUMVAR');
  
  % check for emptyness: cures Matlab bug
  % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
  % throughout this routine use e.g.  sparse(ib,jb,reshape(D(jb),size(sb)).*sb,N,m);
  % rather than  sparse(ib,jb,D(jb).*sb(:),N,m);
  
  if ~isa(a,'gradient')              % non-gradient (scalar or array) ./ gradient array
    if issparse(b.dx)
      m = prod(size(b.x));
      r.x = a ./ b.x;
      D = -r.x ./ b.x;
      [ib,jb,sb] = find(b.dx);
      r.dx = sparse(ib,jb,reshape(D(ib),size(sb)).*sb,m,N);
    else
      r.x = a ./ b.x;
      D = r.x ./ b.x;
      r.dx = b.dx .* repmat(-D(:),1,N);
    end
  elseif ~isa(b,'gradient')          % gradient ./ non-gradient array
    r.x = a.x ./ b;
    bb = 1 ./ b(:);
    if scalar_a                     % gradient scalar ./ non-gradient array
      r.dx = bb * a.dx;
    else                            % gradient array ./ non-gradient array
      if issparse(a.dx)
        m = prod(size(a.x));
        [ia,ja,sa] = find(a.dx);
        r.dx = sparse(ia,ja,reshape(bb(ia),size(sa)).*sa,m,N);
      else
        r.dx = repmat(bb,1,N) .* a.dx;
      end
    end
  else                              % gradient ./ gradient array
    m = prod(size(b.x));
    r.x = a.x ./ b.x;
    if scalar_a                     % gradient scalar ./ gradient array
      if issparse(a.dx) | issparse(b.dx)
        bx = sparse(b.x(:));
        D = 1 ./ sqr(bx);
        Num = bx * a.dx - a.x * b.dx;
        [i,j,s] = find(Num);
        r.dx = sparse(i,j,s.*reshape(D(i),size(s)),m,N);
      else
        Num = b.x(:) * a.dx - a.x * b.dx;
        r.dx = repmat( sqr((1./b.x(:))) , 1 , N ) .* Num;
      end
    else                            % gradient array ./ gradient array
      if issparse(a.dx) | issparse(b.dx)
        r.dx = sparse(1:m,1:m,1./sqr(b.x(:))) * ...   % thanks to J. Kubitz
          ( sparse(1:m,1:m,b.x(:))*a.dx - sparse(1:m,1:m,a.x(:))*b.dx ) ;
      else
        Num = repmat(b.x(:),1,N) .* a.dx - b.dx .* repmat(a.x(:),1,N);
        r.dx = repmat( sqr((1./b.x(:))) , 1 , N ) .* Num;
      end
    end
  end

  r = class(r,'gradient');
  
  if rndold
    setround(rndold)
  end
