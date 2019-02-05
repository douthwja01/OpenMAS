function [e,c] = collect(e,c)
%COLLECT      Collect common exponents in e for multivariate polynomials
%
%On input,
% e   m x n array of exponents, possibly with rows occuring several times
% c   m x 1 array of corresponding (polynomial) coefficients
%
%On output, identical rows of exponents are summed up into corresonding
%  coefficient c
%

% written  08/31/00     S.M. Rump
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  if all(c==0)              % zero polynomial
    e = zeros(1,size(e,2));
    c = typeadj(0,typeof(c));
    return
  end
  
  if ~any(e)                % all zero exponents
    e = e(1,:);
    c = sum(c);
    return
  end
  
  ee = 1e-30;
  if 1+ee==1-ee                         % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  emax = max(e,[],1)+1; 
  index = exp2index(e,emax)+1;
  if any(index>=2^31)           % maximum index for sparse matrices
    indexc_0 = ( c==0 ); 
    if any(indexc_0)
      e(indexc_0,:) = [];
    end
    [e,sortindex] = sortrows(e);
    esortdiff = diff(e,1,1);    
    indexe = logical([1;any(esortdiff,2)]);
    index = cumsum(indexe,1);
    if isa(c,'intval')      
      if any(indexc_0)
        c = intval(c.inf(~indexc_0),c.sup(~indexc_0),'infsup');
      end
      c = intval(c.inf(sortindex),c.sup(sortindex),'infsup');
      i = find( sparse( index , 1 , c~=0 ) );
      setround(-1)
      cinf = sparse( index , 1 , c.inf );
      setround(1)
      csup = sparse( index , 1 , c.sup );
      c = infsup(full(cinf(i)),full(csup(i)));
    else      
      if any(indexc_0)
        c(indexc_0) = [];
      end
      c = c(sortindex);
      [i,j,c] = find( sparse( index , 1 , c ) );      
      if isempty(c)            % zero polynomial
        e = zeros(1,length(emax));
        c = 0;
        setround(rndold)
        return
      end      
    end
    index = 1:size(e,1);
    e = e(index(indexe),:);
  else
    if isa(c,'intval')
      i = find( sparse( index , 1 , double(c~=0) ) );
      setround(-1)
      cinf = sparse( index , 1 , c.inf );
      setround(1)
      csup = sparse( index , 1 , c.sup );
      c = infsup(full(cinf(i)),full(csup(i)));
    else
      [i,j,c] = find( sparse( index , 1 , c ) );
      if isempty(c)            % zero polynomial
        e = zeros(1,length(emax)); 
        c = 0;
        setround(rndold)
        return
      end
    end
    e = index2exp(i-1,emax); 
  end

  setround(rndold)
