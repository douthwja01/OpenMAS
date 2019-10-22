function r = erf(a)
%ERF          Hessian (elementwise) error function
%

% written  31/05/13     S.M. Rump
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  a.x = full(a.x);
  K = prod(size(a.x));
  % factorLB <= 2/sqrt(pi) <= factorUB
  INTLAB_STDFCTS_ERF = getappdata(0,'INTLAB_STDFCTS_ERF');
  factorLB = INTLAB_STDFCTS_ERF.TWO_SQRTPIINF;  % round to nearest
  factorUB = INTLAB_STDFCTS_ERF.TWO_SQRTPISUP;  % ~ 1.12

  if K==1                   % scalar hessian
    
    ax = exp(-a.x(:).^2);
    if isa(a.x,'intval')
        ax = intval(factorLB,factorUB,'infsup') * ax;
    else
      ax = factorLB * ax;
    end
    r.x = erf(a.x);
    r.dx = ax * a.dx;
    r.hx = ax * ( a.hx - reshape( (a.x.*a.dx) * a.dx.' , size(a.hx) ) );
    
  else                      % matrix hessian
    
    N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');
    N2 = N^2;
    
    r.x = erf(a.x);
    a.x = a.x(:);
    if issparse(a.hx)               % input sparse
      
      ax = exp(-a.x.^2);
      if isa(a.x,'intval')
        ax = intval(factorLB,factorUB,'infsup') * ax;
      else
        ax = factorLB * ax;
      end
      sizeax = length(ax);
      [ia,ja,sa] = find(a.dx);
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if isempty(ia)
        r.dx = sparse([],[],[],N,sizeax);
        r.hx = a.hx;
      else
        if isa(a.x,'intval')          % sparse intval
          rdx = times(ax(ja),sa(:),0);
          adx1 = times(a.x(ja),sa(:),0);
          r.dx = intval( sparse(ia,ja,rdx.inf,N,sizeax) , sparse(ia,ja,rdx.sup,N,sizeax) , 'infsup' );
          rdx = intval(sparse(ia,ja,adx1.inf,N,sizeax),sparse(ia,ja,adx1.sup,N,sizeax),'infsup');
        else                          % sparse point  
          r.dx = sparse(ia,ja,ax(ja).*sa(:),N,sizeax); 
          rdx = sparse(ia,ja,a.x(ja).*sa(:),N,sizeax);
        end
        r.hx = a.hx - adx2rhx(N,sizeax,rdx,a.dx);
      end
      [ia,ja,sa] = find(r.hx);        % sparse point or intval
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if ~isempty(ia)
        if isa(ax,'intval')
          rhx = times(ax(ja),sa(:),0);
          r.hx = sparse(ia,ja,intval(rhx.inf,rhx.sup,'infsup'),N2,sizeax);
        else
          r.hx = sparse(ia,ja,ax(ja).*sa(:),N2,sizeax);
        end
      end
      
    else                            % input full
      
      if isa(a.x,'intval')
        ax = intval(factorLB,factorUB,'infsup') * exp((-((a.x).').^2));
      else
        ax = factorLB * exp((-((a.x).').^2));
      end
      ax = ax(ones(N*N,1),:);
      r.dx = a.dx .* ax(1:N,:);
      adx = repmat(a.x(:).',N,1) .* a.dx;
      r.hx = ( a.hx - adx(repmat(1:N,N,1),:) .* a.dx(repmat(1:N,1,N),:) ) .* ax;
      
    end
    
  end
  
  r = class(r,'hessian');
  
  if rndold
    setround(rndold)
  end
