function r = exp(a)
%EXP          Hessian (elementwise) exponential
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

  K = prod(size(a.x));
  if K==1                   % scalar hessian
    
    r.x = exp(a.x);
    r.dx = r.x * a.dx;
    r.hx = r.x * ( a.hx + reshape( (0.5*a.dx) * a.dx.' , size(a.hx) ) );
    
  else                      % matrix hessian
    
    N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');
    N2 = N^2;
    
    r.x = exp(a.x);
    if issparse(a.hx)               % input sparse
      
      ax = r.x(:);
      sizeax = length(ax);
      [ia,ja,sa]=find(a.dx);
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if isempty(ia)
        r.dx = sparse([],[],[],N,sizeax);
        r.hx = sparse([],[],[],N2,sizeax);
      else
        if isa(a.x,'intval')          % sparse intval
          rdx = times(ax(ja),sa(:),0);
          if rdx.complex
            r.dx = intval( sparse(ia,ja,rdx.mid,N,sizeax) , sparse(ia,ja,rdx.rad,N,sizeax) , 'midrad' );
          else
            r.dx = intval( sparse(ia,ja,rdx.inf,N,sizeax) , sparse(ia,ja,rdx.sup,N,sizeax) , 'infsup' );
          end
        else                          % sparse point  
          r.dx = sparse(ia,ja,ax(ja).*sa(:),N,sizeax);        
        end                           
        r.hx = adx2rhx(N,sizeax,0.5*r.dx,a.dx);
      end
      [ia,ja,sa] = find(a.hx);        % sparse point or intval
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if ~isempty(ia)
        if isa(a.x,'intval')
          rhx = times(ax(ja),sa(:),0);
          if rhx.complex
            r.hx = r.hx + intval( sparse(ia,ja,rhx.mid,N2,sizeax) , sparse(ia,ja,rhx.rad,N2,sizeax) , 'midrad' );
          else
            r.hx = r.hx + intval( sparse(ia,ja,rhx.inf,N2,sizeax) , sparse(ia,ja,rhx.sup,N2,sizeax) , 'infsup' );
          end
        else
          r.hx = r.hx + sparse(ia,ja,ax(ja).*sa(:),N2,sizeax);
        end
      end
      
    else                            % input full
      
      r.x = exp(a.x);
      rx = r.x(:).';
      rx = rx(ones(N*N,1),:);
      r.dx = rx(1:N,:) .* a.dx;
      adx = 0.5*a.dx;
      r.hx = rx .* ( a.hx + adx(repmat(1:N,N,1),:) .* a.dx(repmat(1:N,1,N),:) );
      
    end
    
  end
  
  r = class(r,'hessian');
  
  if rndold
    setround(rndold)
  end
