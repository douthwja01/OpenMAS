function res = binom(n,k)
%BINOM        Binomial coefficient, either n or k may be vector, exact for n<=57
%
%   res = binom(n,k)
%
%Result is row vector in case n or k is vector
%

% written  11/05/94     S.M. Rump
% modified 10/13/99     S.M. Rump  vector input for n allowed
% modified 07/16/02     S.M. Rump  vector input k allowed, use of global variable, interval input
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  sizen = size(n);
  sizek = size(k);
  vectorn = ( prod(sizen)~=1 );
  vectork = ( prod(sizek)~=1 );
  if vectorn & vectork
    error('not both arguments can be vector in binom')
  end
  
  IV = 0;
  if isa(n,'intval') 
    IV = 1;
    if ( diam(n)~=0 ) | ( round(inf(n))~=inf(n) )
      res = gamma(n+1)./gamma(k+1)./gamma(n-k+1);  
      setround(rndold)
      return
    end
    n = inf(n);
  end
  if isa(k,'intval') 
    IV = 1;
    if ( diam(k)~=0 ) | ( round(inf(k))~=inf(k) )
      res = gamma(n+1)./gamma(k+1)./gamma(n-k+1);  
      setround(rndold)
      return
    end
    k = inf(k);
  end
  
  if min(n)<max(k)
    error('first argument must be greater equal to second argument in binom')
  end
  
  Nmax = 60;
  INTLAB_BINOM = getappdata(0,'INTLAB_BINOM');	% INTLAB_BINOM(n,k) = (n-1) nchoose (k-1)
  if isempty(INTLAB_BINOM)
    INTLAB_BINOM(4,Nmax) = 0;
    INTLAB_BINOM(1:4,1:4) = [1 0 0 0;1 1 0 0;1 2 1 0;1 3 3 1];
    setappdata(0,'INTLAB_BINOM',INTLAB_BINOM);  % update INTLAB_BINOM
  end
  
  nmax = max(n);
  Gmax = size(INTLAB_BINOM,1);
  if Gmax<(nmax+1)
    for i=(Gmax+1):min(nmax+1,Nmax+1)
      v = INTLAB_BINOM(i-1,1:i-1);
      INTLAB_BINOM(i,1:i) = [0 v] + [v 0];
    end
    setappdata(0,'INTLAB_BINOM',INTLAB_BINOM);  % update INTLAB_BINOM
  end
  
  if prod(sizen)>1            % n is vector
    if k>Nmax                 % k is big
      res = bigbinom(n,k,vectorn,vectork,IV);
    else                      % k is small
      index = ( n>Nmax );
      if any(index)           % some or all n big        
        res = INTLAB_BINOM((min(n)+1):(Nmax+1),k+1);
        res = [res(:)' bigbinom(n(index),k,vectorn,vectork,IV)];
      else
        res = INTLAB_BINOM(n+1,k+1);
        res = res(:)';
      end
    end
  else                        % n is scalar, k may be vector
    if n>Nmax                 % n is big
      res = bigbinom(n,k,vectorn,vectork,IV);
    else                      % n is small, k must be small
      index = ( k>Nmax );
      res = INTLAB_BINOM(n+1,k+1);
      res = res(:)';
    end
  end

  if vectorn
    res = reshape(res,sizen);
  end
  
  if vectork
    res = reshape(res,sizek);
  end
    
  setround(rndold)
 
  
  function res = bigbinom(n,k,vectorn,vectork,IV)
  %internal for big values of n, res = gamma(n+1)./gamma(k+1)./gamma(n-k+1);
  
  if IV                       % interval computation
    if vectorn                % n is vector, k must be scalar
      Num = ones(k,length(n));
      Num(1,:) = n-k+1;
      setround(-1)
      Numinf = cumprod(cumsum(Num));
      Numinf = Numinf(k,:);
      setround(1)
      Numsup = cumprod(cumsum(Num));
      Numsup = Numsup(k,:);
      Densup = prod(1:k);
      setround(-1)
      resinf = Numinf./Densup;
      Deninf = prod(1:k);
      ressup = Numsup./Deninf;
    elseif vectork            % k is vector, n must be scalar
      setround(-1)
      Numinf = cumprod(n:-1:(n-max(k)+1));
      Deninf = cumprod(1:max(k));
      setround(1)
      Numsup = cumprod(n:-1:(n-max(k)+1));
      Densup = cumprod(1:max(k));
      ressup = Numsup(k)./Deninf(k);
      setround(-1)
      resinf = Numinf(k)./Densup(k);
    else                      % neither n nor k vector
      setround(-1)
      Numinf = prod((n-k+1):n);
      Deninf = prod(1:k);
      setround(1)
      Numsup = prod((n-k+1):n);
      Densup = prod(1:k);
      ressup = Numsup./Deninf;
      setround(-1)
      resinf = Numinf./Densup;
    end
    res = intval(resinf,ressup,'infsup');
  else
    if vectorn                % n is vector, k must be scalar
      Num = ones(k,length(n));
      Num(1,:) = n-k+1;
      Num = cumprod(cumsum(Num,1),1);
      Num = Num(k,:);
      res = Num./prod(1:k);
    elseif vectork            % k is vector, n must be scalar
      Num = cumprod(n:-1:(n-max(k)+1));
      Den = cumprod(1:max(k));
      res = Num(k)./Den(k);
    else                      % neither n nor k vector
      Num = prod((n-k+1):n);
      res = Num./prod(1:k);
    end
  end
