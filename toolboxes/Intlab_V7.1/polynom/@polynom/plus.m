function r = plus(p,q)
%PLUS         Polynomial addition  p + q
%
%p or q may be scalar (interval)
%
%Note that the result is always of type 'polynom', even if the result is constant.
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(p,'polynom')        % p is scalar
    if p==0
      r = q;  
      if rndold
        setround(rndold)
      end
      return
    end
    if size(q.e,2)<=1         % q univariate
      r = q;
      m = length(q.c);
      r.c = typeadj( r.c , typeof(p) );
      r.c(m) = r.c(m) + p;
    else                      % q multivariate
      r = q;
      index = find( sum(q.e,2)==0 );
      if isempty(index)       % constant term not in polynomial q
        r.e = [ zeros(1,size(r.e,2)) ; r.e ];
        r.c = typeadj( [ 0 ; r.c ] , typeof(p) );
        r.c(1) = p;
      else
        r.c = typeadj( r.c , typeof(p) );
        r.c(index) = r.c(index) + p;
        if iszero(r.c(index))
          r.e(index,:) = [];
          r.c(index) = [];
        end
      end
    end
    r = normalize(r);  
    if rndold
      setround(rndold)
    end
    return
  end

  if ~isa(q,'polynom')        % q is scalar
    if q==0
      r = p;  
      if rndold
        setround(rndold)
      end
      return
    end
    if size(p.e,2)<=1         % p univariate
      r = p;
      m = length(p.c);
      r.c = typeadj( r.c , typeof(q) );
      r.c(m) = r.c(m) + q;
    else                      % p multivariate
      r = p;
      index = find( sum(p.e,2)==0 );
      if isempty(index)       % constant term not in polynomial p
        r.e = [ zeros(1,size(r.e,2)) ; r.e ];
        r.c = typeadj( [ 0 ; r.c ] , typeof(q) );
        r.c(1) = q;
      else
        r.c = typeadj( r.c , typeof(q) );
        r.c(index) = r.c(index) + q;
        if iszero(r.c(index))
          r.e(index,:) = [];
          r.c(index) = [];
        end
      end
    end
    r = normalize(r);  
    if rndold
      setround(rndold)
    end
    return
  end
  
  if isempty(p.v)                % p constant
    p.e = zeros(1,size(q.e,2));
    p.v = q.v;
  elseif isempty(q.v)            % q constant
    q.e = zeros(1,size(p.e,2));
    q.v = p.v;
  end

  if ( size(p.e,2)<=1 ) & ...       % p univariate or constant
     ( size(q.e,2)<=1 ) & ...       % q univariate or constant
       isequal(p.v,q.v)             % variable of p and q identical
    np = length(p.c);
    nq = length(q.c);
    r.e = max(np,nq)-1;
    if np<nq
      r.c = [zeros(1,nq-np) p.c] + q.c;
    elseif nq<np
      r.c = p.c + [zeros(1,np-nq) q.c];
    else
      r.c = p.c + q.c;
    end;
    r.v = p.v;
    r = class(r,'polynom');
    r = normalize(r);  
    if rndold
      setround(rndold)
    end
    return
  end

  % result p+q multivariate
  if isequal(p.v,q.v)                   % identical set of variables
    r.e = [ p.e ; q.e ];                % both p and q must be multivariate
    r.c = [ p.c ; q.c ];
    r.v = p.v;
  else
    np = length(p.c);
    nq = length(q.c);
    [z,I,J] = joinvars(p.v,q.v); 
    r.e = zeros(np+nq,length(z));
    if size(p.e,2)>1                    % p multivariate
      r.e(1:np,I) = p.e;
    else                                % p univariate
      r.e(1:np,I) = (np-1:-1:0)';
    end
    if size(q.e,2)>1                    % q multivariate
      r.e(np+1:end,J) = q.e;
    else                                % q univariate
      r.e(np+1:end,J) = (nq-1:-1:0)';
    end
    r.c = [ p.c(:) ; q.c(:) ];          % (:) covers univariate polynomials
    r.v = z;
  end

  r = normalize( class(r,'polynom') );
  
  if rndold
    setround(rndold)
  end
