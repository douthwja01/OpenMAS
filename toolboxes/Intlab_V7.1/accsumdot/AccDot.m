function [res,exact] = AccDot(varargin)
%AccDot       Accurate dot product with faithful, to nearest or K-fold rounding
%
%   res = AccDot(A1,B1,...,An,Bn,K)
%
%res is sum(Ai*Bi)
%   K optional, default K=1 for faithfully rounded result
%   K   =0    rounded to nearest result
%       >1    result in cell array of length K of non-overlapping sequences,
%               sum(res{i})=sum(Ai,Bi) with K-fold accuracy
%       =inf  result in cell array of length K of non-overlapping sequences of
%               sufficient length s.t. sum(res{i})=sum(Ai,Bi) exactly
%       []    interval inclusion of the result; inclusion is best possible if no
%               underflow occurs.
%
%output exact (optional), array of 0/1 of size of result with entry equal to 1
%  iff corresponding dot product exact.
%
%All factors may be real non-interval. 
%  The size of every product Ai*Bi must be the same.
%
%One of the factors Ai,Bi may be a cell array. If Ai is a cell array, say,
%  then Ai*Bi is interpreted as sum(Ai{k})*Bi. All entries Ai{k} must be of 
%  the same size. This gives a convenient way to compute an approximate inverse
%  of an extremely ill-conditioned matrix only using double precision and our
%  K-fold accurate dot product AccDot. For details and examples with cond(A)~1e65
%  and more, see the INTLAB homepage.
%
%Interval input makes only sense, if all intervals are degenerated (point intervals); 
%  therefore omitted.
%
%All factors are first stored into one array to be summed up; summation of 
%  components by loops. For matrix input Ai,Bi, there may be interpretation overhead.
%For speed and to fight interpretation overhead, 2 extra matrices are necessary; for
%  the product of to nxn matrices this would be n^3 memory, therefore this is not allowed.
%
%   [res,exact] = AccDot(A1,B1,...,An,Bn,K)
%
%output exact (optional), array of 0/1 of size of result with entry equal to 1
%  iff corresponding dot product exact.
%
%Result is set to NaN if overflow occurs. For simplicity, factors Ai,Bi and 
%  products Ai*Bi are limited in absolute value to less or equal 1e300. Small
%  numbers and underflow is treated correctly, rounding to nearest might be off
%  one bit if underflow occurs.
%
%A simple way to compute an inclusion of sum(Ai*Bi) of K-fold accuracy is
%
%   [res,exact] = AccDot(A1,B1,...,An,Bn,K);
%   ResK = intval(res{K});
%   ResK(~exact) = infsup(pred(ResK(~exact)),succ(ResK(~exact)));
%   res{K} = ResK;
%
%This uses faithfully rounded result. Note that the interval is represented by
%   sum(res(1:K-1)) + ResK;
%This interval is generally 2 bits wide unless res(mu,nu)=sum(Ai*Bi)(mu,nu). 
%
%A best possible interval including sum(Ai*Bi) is computed by
%
%   res = AccDot(A1,B1,...,An,Bn,[]);
%
%The inclusion is best possible, however, more expensive than the previous approach. 
%The gain is at most one bit.
%Maximum number of nonzero elements per sum is limited to 67108862, which
%seems sufficient for Matlab. More elements can be treated, see our paper (Huge length).
%
%Uses various algorithms in
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, SIAM J. Sci. Comput., 31(1):189-224, 2008.
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%
%CAUTION: !!! THIS IMPLEMENTATION SUFFERS SEVERELY FROM INTERPRETATION OVERHEAD !!!
%!!! IT IS INCLUDED TO SHOW THE PRINCIPLES OF THE NEW METHODS !!!
%!!! DO NOT USE FOR LARGE DIMENSIONS !!!
%

% written  12/12/05     S.M. Rump
% modified 02/13/06     S.M. Rump  check for complex operands
% modified 08/26/12     S.M. Rump  rounding
% modified 12/23/12     S.M. Rump  typo
%

  e = 1e-30;
  if 1+e==1-e                       % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  % detect parameters
  len = length(varargin);
  if even(len)                          % default accuracy
    K = 1;
  else
    K = varargin{len};
  end
  
  % check K
  if ~(K>=0) & ~isempty(K)
    error('invalid value K for accuracy in AccDot')
  end
    
  % detect dimension of result
  Ai = varargin{1};
  if iscell(Ai)
    Ai = Ai{1};
  end
  Bi = varargin{2};
  if iscell(Bi)
    Bi = Bi{1};
  end
  if prod(size(Ai))==1
    [m n] = size(Bi);
  elseif prod(size(Bi))==1
    [m n] = size(Ai);
  else
    m = size(Ai,1);                     % m or n must be 1
    n = size(Bi,2);
  end
  sizeres = [m n];
    
  % Convert input into A*B+C, where A*B and C is m x n  
  % C covers products Ai*Bi where one factor is +/-1, stored as index set indexC
  % index is negative if cofactor is cell array
  A = [];
  B = [];
  indexC = [];                          % summands, cofactor +/- 1
  indexD = [];                          % cell array summands, cofactor +/- 1
  isfull = 1;                           % all factors full
  
  for i=1:floor(nargin/2)    
    % next product Ai*Bi
    Ai = varargin{2*i-1};
    Bi = varargin{2*i};
    if iscell(Ai) & iscell(Bi)
      error('at most factor can be cell array in AccDot')
    end
    if iscell(Ai)      
      isfull = isfull & ( ~issparse(Ai{1}) & ~issparse(Bi) );
      sAi = size(Ai{1});
      for k=2:length(Ai)
        if ~isequal(sAi,size(Ai{k}))
          error('sizes in cell array factors in AccDot must coincide')
        end
      end
      sBi = size(Bi);
    elseif iscell(Bi)      
      isfull = isfull & ( ~issparse(Ai) & ~issparse(Bi{1}) );
      sAi = size(Ai);
      sBi = size(Bi{1});
      for k=2:length(Bi)
        if ~isequal(sBi,size(Bi{k}))
          error('sizes in cell array factors in AccDot must coincide')
        end
      end
    else   
      isfull = isfull & ( ~issparse(Ai) & ~issparse(Bi) );
      sAi = size(Ai);
      sBi = size(Bi);
    end
    % check dimensions
    if prod(sAi)==1 
      if ~isequal(sBi,sizeres)
      error('dimensions in AccDot do not match')
      end
    elseif  prod(sBi)==1
      if ~isequal(sAi,sizeres)   
        error('dimensions in AccDot do not match')
      end
    elseif sAi(2)~=sBi(1)
      error('dimensions in AccDot do not match')
    end
    
    % sAi size of factor Ai,  sBi size of factor Bi
    % determine factors
    if iscell(Ai)                         % first factor Ai cell array
      if prod(sBi)==1                     % second factor Bi scalar
        if abs(Bi)==1                     % Bi = +/- 1
          indexD = [ indexD Bi*(2*i-1) ];
        else                              % Bi scalar, but not +/- 1
          if isfull
            Bi_ = sparse(1:n,1:n,Bi);
          else
            Bi_ = Bi*eye(n);
          end
          for k=1:length(Ai)
            A = [ A Ai{k} ];
            B = [ B ; Bi_ ];
          end
        end
      else                                % second factor Bi non-scalar
        if ( prod(sAi)==1 ) & ( ~isequal(sizeres,[1 1]) )
          error('cell array factor cannot be sequence of scalars')
        end
        for k=1:length(Ai)
          A = [ A Ai{k} ];
          B = [ B ; Bi ];
        end
      end
    elseif iscell(Bi)                     % second factor Bi cell array
      if prod(sAi)==1                     % first factor Ai scalar
        if abs(Ai)==1                     % Ai = +/- 1
          indexD = [ indexD Ai*(2*i) ];
        else                              % Ai scalar, but not +/- 1
          if isfull
            Ai_ = sparse(1:m,1:m,Ai);
          else
            Ai_ = Ai*eye(m);
          end
          for k=1:length(Bi)
            A = [ A Ai_ ];
            B = [ B ; Bi{k} ];
          end
        end
      else                                % first factor Ai non-scalar
        if ( prod(sBi)==1 ) & ( ~isequal(sizeres,[1 1]) )
          error('cell array factor cannot be sequence of scalars')
        end
        for k=1:length(Bi)
          A = [ A Ai ];
          B = [ B ; Bi{k} ];
        end
      end
    else                                  % none of Ai,Bi cell array
      if prod(sAi)==1                     % first factor Ai scalar
        if abs(Ai)==1                     % Ai = +/- 1
          indexC = [ indexC Ai*(2*i) ];
        else
          if isfull
            A = [ A Ai*eye(m) ];
          else
            A = [ A sparse(1:m,1:m,Ai) ]; % Ai scalar, but not +/- 1
          end
          B = [ B ; Bi ];
        end
      elseif prod(sBi)==1                 % second factor Bi scalar
        if abs(Bi)==1                     % Bi = +/- 1
          indexC = [ indexC Bi*(2*i-1) ];
        else
          A = [ A Ai ];                   % Bi scalar, but not +/- 1
          if isfull
            B = [ B ; sparse(1:n,1:n,Bi) ];
          else
            B = [ B ; Bi*eye(n) ];
          end
        end
      else                                % both factors Ai,Bi non-scalar
        A = [ A Ai ];
        B = [ B ; Bi ];
      end
    end
  end
  
  % check interval input
  if isa(A,'intval') | isa(B,'intval')
    error('input of AccDot must not be of type intval')
  end
  
  % check for complex input
  if ( ~isreal(A) ) | ( ~isreal(B) )
    error('input of AccDot must be real')
  end
      
  % check huge input
  if any(any(abs(A)>1e300)) | any(any(abs(B)>1e300))
    error('for simplicity all entries and intermediate products limited to 1e300')
  end
  if ~isempty(indexC)
    for k=1:length(indexC)      
      if any(any(varargin{abs(indexC(k))}>1e300))
        error('for simplicity all entries and intermediate products limited to 1e300')
      end
    end
  end
  kD = 0;
  if ~isempty(indexD)
    for k=1:length(indexD)
      DD = varargin{abs(indexD(k))};
      kD = kD + length(DD);
      for kk=1:length(DD)
        if any(any(DD{kk}>1e300))
          error('for simplicity all entries and intermediate products limited to 1e300')
        end
      end
    end
  end

  % transpose A:  result A.'*B
  A = A.';
  
  factor = 134217729;               % splitting factor 2^27+1
  % [A1,A2] = Split(A)
  C = factor*A;
  A1 = C - ( C - A );               % upper part of A
  A2 = A - A1;                      % A = A1+A2 exact splitting
  % [B1,B2] = Split(B)
  C = factor*B;
  B1 = C - ( C - B );               % upper part of B
  B2 = B - B1;                      % B = B1+B2 exact splitting
  
  % initialize result
  if isempty(K)
    rescell = 0;
  else
    rescell = ( K>1 );
  end
  if isfull
    if rescell
      res{1} = zeros(m,n);
    else
      res = zeros(m,n);
    end
  else
    if rescell
      res{1} = sparse([],[],[],m,n);
    else
      res = sparse([],[],[],m,n);
    end
  end
  if ( K>1 ) & ~isinf(K)
    for k=2:K
      res{k} = res{1};
    end
  end
  if isempty(K)
    res = intval(res);
  end
  if nargout==2
    if isfull
      exact = zeros(m,n);
    else
      exact = sparse([],[],[],m,n);
    end
  end
  
  % compute result
  sA = size(A,1);
  sC = length(indexC);
  sD = length(indexD);
  summands = ( sC+sD~=0 );
  % underflow constants
  phi = 2^643;                                % factor sqrt(eps^-4 eta^-1) for Ai,Bi
  Fbound = 2^(-968);                          % lower bound eps^-2 eta for factors
  for i=1:m
    for j=1:n
      if sA~=0
        x1 = A1(:,i);
        x2 = A2(:,i);
        y1 = B1(:,j);
        y2 = B2(:,j);
        AiBip = A(:,i).*B(:,j);
        r = x2.*y2 - ((( AiBip - x1.*y1 ) - x2.*y1 ) - x1.*y2 );
        if any(abs(AiBip)>1e300)              % AiBip is vector
          error('for simplicity all entries and intermediate products limited to 1e300')
        end
      else
        AiBip = [];
        r = [];
      end
      if summands
        Ci = zeros(sC+kD,1);
        for k=1:sC
          if indexC(k)>0
            Ci(k) = varargin{indexC(k)}(i,j);
          else
            Ci(k) = -varargin{-indexC(k)}(i,j);
          end
        end
        kkk = 0;
        for k=1:sD
          DD = varargin{abs(indexD(k))};
          if indexD(k)>0
            for kk=1:length(DD)
              kkk = kkk+1;
              Ci(sC+kkk) = DD{kk}(i,j);
            end
          else
            for kk=1:length(DD)
              kkk = kkk+1;
              Ci(sC+kkk) = -DD{kk}(i,j);
            end
          end
        end
      else
        Ci = [];
      end
      p = [ AiBip ; Ci ];
      % execute one ExtractVector
      np = nnz(p);                          % initialization
      mu = full(max(abs(p)));               % abs(p_i) <= mu; full: avoid matlab bug
      if ( np==0 ) | ( mu==0 )              % no or only zero summands
        R = 0;                              % result exactly zero in any rounding
        p = [];
      else
        Ms = 2^nextpow2(np+2);              % np+2 <= 2^M
        if Ms^2*eps>1
          error('vector length n too large for AccDot; Huge n not implemented')
        end
        sigma = Ms*2^nextpow2(mu);          % first extraction unit
        if isinf(sigma) | isnan(sigma)      % overflow, could be avoided with scaling
          if iscell(res)
            res{1}(i,j) = NaN;
          else
            res(i,j) = NaN;
          end
          exact(i,j) = 0;
          continue
        else
          q = ( sigma + p ) - sigma;        % [R,p] = ExtractVector(sigma,p);
          R = sum(q);                       % sum of leading terms
          p = p - q;                        % remaining terms
        end
      end
      if K<=1                               % nearest of faithful rounding
        [res(i,j),ext,R,p] = AccSum( [ p ; r ] , K , 0 , R );
        if res(i,j)==realmin                % take care of underflow, original sum = R + sum(p)
          index = find( abs(AiBip)<Fbound );      % small parts of AiBip
          p = [ p ; -AiBip(index) ; -r(index) ];  % subtract small parts
          [res(i,j),ext] = ...
            AccSum( scaleinput(A(index,i),B(index,j),p,phi) , K , 0 , phi*(phi*R) );
          if ( nargout==2 ) & ext
            exact(i,j) = diam( (intval(res(i,j))/phi) / phi == 0 );
          end
          res(i,j) = (res(i,j)/phi) / phi;  % avoid overflow
        end
      elseif isempty(K)                       % interval result
        [rr,ext,R,p] = AccSum( [ p ; r ] , 0 , -1 , R );
        if rr==realmin                        % take care of underflow, original sum = R + sum(p)
          index = find( abs(AiBip)<Fbound );  % small parts of AiBip
          p = [ p ; -AiBip(index) ; -r(index) ];  % subtract small parts
          [rr,ext] = AccSum( scaleinput(A(index,i),B(index,j),p,phi) , 0 , -1 , phi*(phi*R) );
          if ext
            res(i,j) = intval(rr);
          else
            res(i,j) = infsup(rr,succ(rr));
          end
          res(i,j) = ( res(i,j)/phi ) / phi;
        else
          res(i,j) = infsup(rr,succ(rr));
        end
      else                                          % K-fold accuracy
        p = [ p ; r ];
        scale = 0;
        cont = 1;
        k = 0;
        ext(i,j) = 0;
        while ( k<K ) & cont
          k = k+1;
          [rr,ext,R,p] = AccSum( p , 1 , 0 , R );
          if ~scale & ( rr==realmin )               % take care of underflow, original sum = R + sum(p)
            scale = 1;                              % only one scaling
            index = find( abs(AiBip)<Fbound );      % small parts of AiBip
            p = [ p ; -AiBip(index) ; -r(index) ];  % subtract small parts
            [rr,ext,R,p] = ...                      % scaled result
              AccSum( scaleinput(A(index,i),B(index,j),p,phi) , 1 , 0 , phi*(phi*R) );
          end
          rrr = rr;
          if scale
            rr = ( rr/phi )/phi;
          end
          if isinf(K) & ( k>length(res) ) & ( rr~=0 )
            if isfull
              res{k} = zeros(m,n);
            else
              res{k} = sparse([],[],[],m,n);
            end
          end
          if rr~=0
            res{k}(i,j) = rr;
          end
          if abs(rr)<=realmin                       % result exact
            cont = 0;
          end
          if isnan(rr)
            res{k}(i,j) = NaN;
            cont = 0;
          end
          if ( nargout==2 ) & ext
            exact(i,j) = ( diam( (intval(rrr)/phi) / phi )==0 );
          end
        end
      end
    end
  end

  if rndold
    setround(rndold)
  end

  
function [res,exact,R,p] = AccSum(p,K,rnd,rho)
%Simplified version of AccSum :
%  rnd only valid for K==0
%  no overflow can occur
%  indices for possible underflow stored
%  vector length not too large
%  stops if sigma too small
%  ~(K>1)
%

  res = 0;
  exact = 1;
  R = 0;
  if isempty(p)
    return
  end
  
  % input real, compute sum(p)  
  % rounding to nearest
  if K==0                             % rnd maybe -1 or 0 or +1
    % let S:=sum(p)
    [res,dummy,R,p] = AccSum(p,1,0,rho);   % S = res + R + sum(p)  for new p
    if res==realmin                   % take care of underflow
      return
    end
    [delta,dummy,R,p] = AccSum(p,1,0,R);   % delta + R + sum(p) = S - res  for new R,p
    if delta==realmin                 % scale; otherwise sign(delta) = sign(S-res)
      return
    end
    if delta==0                       % result exact: res=S in any rounding
      return
    end
    exact = 0;                        % result not exact
    
    if rnd~=0                         % rounding downwards or upwards
      if sign(delta)==rnd             % res on wrong side
        if rnd==1
          res = succ(res);
        else
          res = pred(res);
        end
      end
      return
    end
    
    % sign(delta) = sign(S-res)
    % compute nearest neighbor resp in direction sign(delta), rnd=0
    if delta>0
      resp = succ(res);
    else
      resp = pred(res);
    end
    mu = (resp-res)/2;
    if abs(delta)>abs(mu)
      res = resp;
    elseif abs(delta)==abs(mu)
      delta = AccSum(p,1,0,R);
      if delta==0
        res = res + mu;
      elseif sign(delta)==sign(mu)
        res = resp;
      end
    end
    return
  end
  
  % the standard case: real vector input, K=1
  if issparse(p)
    n = nnz(p);                         % initialization
  else
    n = length(p);
  end
  tau1 = 0;
  tau2 = 0;  
  mu = full(max(abs(p)));               % abs(p_i) <= mu
  if ( n==0 ) | ( mu==0 )               % no or only zero summands
    res = rho;                          % result exactly zero in any rounding
    return
  end
  Ms = 2^nextpow2(n+2);                 % n+2 <= 2^M
  sigma = Ms*2^nextpow2(mu);            % first extraction unit
  phi = 2^(-53)*Ms;                     % factor to decrease sigma
  factor = 2*phi*Ms;                    % factor for sigma check

  % underflow constant
  sigmamin = 2^(-968);                  % lower bound eps^-2 eta for sigma
  t = rho;
  while 1
    if sigma<=sigmamin                  % sigma too small, but t+tau exact
      res = realmin;                    % indicates underflow
      R = t;                            % original sum = R + sum(p)
      return
    end
    q = ( sigma + p ) - sigma;          % [tau,p] = ExtractVector(sigma,p);
    tau = sum(q);                       % sum of leading terms
    p = p - q;                          % remaining terms
    tau1 = t + tau;                     % new approximation
    if abs(tau1)>=factor*sigma 
      tau2 = tau - ( tau1 - t );        % [tau1,tau2] = FastTwoSum(t,tau)
      res = tau1 + ( tau2 + sum(p) );   % faithfully rounded final result
      R = tau2 - ( res - tau1 );        % only for K-fold result
      if nargout>=2
        exact = ( R==0 ) & ~any(p);
      end
      return
    end
    t = tau1;                           % sum t+tau exact
    if t==0                             % accelerate case sum(p)=0
      [res,exact,R,p] = AccSum(p(p~=0),K,0,0);  % recursive call, zeros eliminated
      return
    end
    sigma = phi*sigma;                  % new extraction unit
  end
    
  
  
function p = scaleinput(Ai,Bi,Ci,phi)
% scale  p = phi*Ai * phi*Bi  for column vector input
 
  factor = 134217729;               % splitting factor 2^27+1
  % [A1,A2] = Split(phi*Ai)
  Ai = phi*Ai;
  C = factor*Ai;
  A1 = C - ( C - Ai );            % upper part of Ai
  A2 = Ai - A1;                   % Ai = A1+A2 exact splitting
  % [B1,B2] = Split(phi*Bi)
  Bi = phi*Bi;
  C = factor*Bi;
  B1 = C - ( C - Bi );            % upper part of Bi
  B2 = Bi - B1;                   % Bi = B1+B2 exact splitting
  % phi^2*A(:,i)*B(:,j) = p+r, no underflow possible
  p = Ai.*Bi;
  p = [ p ; A2.*B2 - ((( p - A1.*B1 ) - A2.*B1 ) - A1.*B2 ) ; phi*(phi*Ci) ];
