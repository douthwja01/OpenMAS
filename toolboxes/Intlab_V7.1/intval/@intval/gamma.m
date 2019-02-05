function y = gamma(x)
%GAMMA        Implements  gamma(x)  for real intervals
%
%   y = gamma(x)
%
%interval standard function implementation
%

% written  06/19/13     S.M. Rump
%

  if x.complex
    error('Gamma function only for real arguments')
  end
  
  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      y = intval(NaN(size(x)));
      index = ~index;
      %VVVV  y(index) = gamma(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,gamma(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = gamma(full(x));
    end
    return
  end
  
  index =  ( ( x.inf<=0 ) & ( ceil(x.inf)<=floor(x.sup) ) ) | ...
    isnan(x.inf) | isnan(x.sup) | ( x.inf==-inf ) | ( x.sup==-inf );
  if any(index)                             % non-negative integers
    y = x;
    y.inf(index) = NaN;
    y.sup(index) = NaN;
    index = ~index;
    yy = gamma(intval(x.inf(index),x.sup(index),'infsup'));
    y.inf(index) = yy.inf;
    y.sup(index) = yy.sup;
    return
  end

  e = 1e-30;
  if 1+e==1-e                               % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  y = x;                                    % x positive or negative
  
  if all( x.inf(:)==x.sup(:) );             % thin input
    
    xx = x.inf(:);
    index = ( xx>=1 );
    if any(index)                           % direct formula for x>=1
      [y.inf(index),y.sup(index)] = gamma_rnd(xx(index),[]);    % lower and upper bound
    end
    
    index = ~index;
    if any(index)                           % treat indices with x<1
      xx = xx(index);
      len = length(xx);                     % number of elements
      col = 1 - floor(xx);
      maxcol = max(col);
      factor = repmat(xx,1,maxcol) + repmat(0:(maxcol-1),len,1);
      factor(factor>1) = 1;
      F = prod(intval(factor),2);
      xx = intval(xx) + col;                % x > 1
      [yinf,ysup] = gamma_rnd(xx.inf,[]);
      index0 = ( xx.inf~=xx.sup );
      if any(index0)
        [y2inf,y2sup] = gamma_rnd(xx.sup(index0),[]);
        yinf(index0) = min(yinf(index0),min(y2inf,y2sup));
        ysup(index0) = max(ysup(index0),max(y2inf,y2sup));
      end
      yy = intval(yinf,ysup,'infsup')./F;
      y.inf(index) = yy.inf;
      y.sup(index) = yy.sup;
    end
    
  else                                      % thick input
    
    xinf = x.inf(:);                        % no non-positive integer in x
    xsup = x.sup(:);
    xmin = hex2num('3ff762d86356be39');     % gamma'(xmin) < 0 < gamma'(succ(xmin))
    gammamin = hex2num('3fec56dc82a74aee'); % gamma(x)>gammamin for x in [1,2]
    
    index = ( xinf<1 );
    if any(index)                           % x < 1
      y1 = gamma(intval(xinf(index)));
      y2 = gamma(intval(xsup(index)));
      y.inf(index) = min(y1.inf,y2.inf);
      y.sup(index) = max(y1.sup,y2.sup);
    end
    
    index = ~index;            
    if any(index)                           % x > 0
      index1 = index & ( xsup<=xmin);       % x left of minimum
      if any(index1)
        y.inf(index1) = gamma_rnd(xsup(index1),-1);
        y.sup(index1) = gamma_rnd(xinf(index1),1);
      end
      index1 = index & ( xmin<xinf);       % x right of minimum
      if any(index1)
        y.inf(index1) = gamma_rnd(xinf(index1),-1);
        y.sup(index1) = gamma_rnd(xsup(index1),1);
      end
      index1 = index & ( xinf<=xmin) & ( xmin<xsup);  % minimum in x
      if any(index1)
        y.inf(index1) = gammamin;
        y.sup(index1) = max(gamma_rnd(xinf(index1),1),gamma_rnd(xsup(index1),1));
      end
    end
        
  end

  y.inf = reshape(y.inf,size(x.inf));
  y.sup = reshape(y.sup,size(x.sup));
  
  setround(rndold)  
