function y = erf(x)
%ERF          Implements  erf(x)  for real intervals
%
%   y = erf(x)
%
%interval standard function implementation
%

% written  05/31/13     S.M. Rump
%

  if x.complex
    error('Error function only for real arguments')
  end
  
  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,erf(full(sx)),m,n);
    return
  end

  e = 1e-30;
  if 1+e==1-e                               % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  y = x;
  
  if all( x.inf(:)==x.sup(:) );             % thin input
    
    xx = x.inf(:);
    index = ( xx>=0 );
    if any(index)                           % direct formula for x>=0
      [y.inf(index),y.sup(index)] = erf_rnd(xx(index),[]);    % lower and upper bound
    end
    
    index = ~index;
    if any(index)
      [yinf,ysup] = erf_rnd(-xx(index),[]);    % lower and upper bound
      dummy = yinf;
      y.inf(index) = -ysup;
      y.sup(index) = -dummy;
    end
    
  else                                      % thick input
    
    xinf = x.inf(:);
    xsup = x.sup(:);
    
    index = ( xinf<=0 );
    if any(index)
      y.inf(index) = - erf_rnd(-xinf(index),1);
    end
    index = ~index; 
    if any(index)
      y.inf(index) = erf_rnd(xinf(index),-1);
    end
    
    index = ( xsup<=0 );
    if any(index)
      y.sup(index) = - erf_rnd(-xsup(index),-1);
    end
    index = ~index; 
    if any(index)
      y.sup(index) = erf_rnd(xsup(index),1);
    end
    
  end

  y.inf = reshape(y.inf,size(x.inf));
  y.sup = reshape(y.sup,size(x.sup));
  
  setround(rndold)  
