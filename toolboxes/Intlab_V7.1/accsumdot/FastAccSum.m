function res = FastAccSum(p)
%FastAccSum   Ultimately fast and accurate summation with faithful rounding
%
%   res = FastAccSum(p)
%
%For real or complex input vector, dense or sparse, the result res is
%sum(p_i) faithfully rounded. Input vector p must not be of type intval.
%
%Maximum number of nonzero elements per sum is limited to 67,108,862 in 
%double precision, which seems sufficient for Matlab.
%
%Implements new algorithm in
%  S.M. Rump: Ultimately Fast Accurate Summation, SIAM Journal on 
%     Scientific Computing (SISC), 31(5):3466-3502, 2009.
%
%CAUTION: !!! THIS IMPLEMENTATION SUFFERS SEVERELY FROM INTERPRETATION OVERHEAD !!!
%!!! IT IS INCLUDED TO SHOW THE PRINCIPLES OF THE NEW METHODS !!!
%!!! DO NOT USE FOR LARGE DIMENSIONS !!!
%

% written  08/28/08     S.M. Rump
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 11/12/12     S.M. Rump  constant (thanks to Heng Tian, Hong Kong)
% modified 12/23/12     S.M. Rump  typo

%

  res = 0;
  if isempty(p)
    return
  end

  % check size
  if length(size(p))>2
    error('FastAccSum not defined for multi-dimensional arrays.')
  end
  p = p(:).';                            % form row vector
  if size(p,1)~=1
    error('FastAccSum only for vector input')
  end

  % check interval input
  if isa(p,'intval')
    error('FastAccSum not defined for interval input')
  end
  
  % check improper input
  if any(isnan(p)) | any(isinf(p))
    res = NaN;
    return
  end

  % take care of complex input
  if ~isreal(p)
    resreal = FastAccSum(real(p));
    resimag = FastAccSum(imag(p));
    res = complex(resreal,resimag);
    return
  end

  % input real, compute sum(p)
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if issparse(p)
    n = nnz(p);                         % initialization
    p = nonzeros(p).';
  else
    n = length(p);
  end
  
  % initialize constants depending on precision
  if isa(p,'single')
    eps = 2^(-24);
    eta = 2^(-149);
  else
    eps = 2^(-53);
    eta = 2^(-1074);
  end
  
  % check dimension
  if ((2*n+4)*n+6)*eps>1
    error('dimension too large for FastAccSum')
  end
  
  % initialize constants
  c1 = 1-n*eps;
  c2 = 1-(3*n+1)*eps;
  c3 = 2*eps;
  c4 = 1-eps;
  c5 = 2*n*(n+2)*eps;
  c6 = 1-5*eps;
  c7 = (1.5+4*eps)*(n*eps);
  c8 = 2*n*eps;
  c9 = eta/eps;
  
  T = sum(abs(p))/c1;                   % sum(abs(p)) <= T
  if T<=c9                              % no rounding error
    res = sum(p);
    if rndold, setround(rndold), end
    return
  end
  tp = 0;
  while 1
    sigma0 = (2*T)/c2;
    P = cumsum([sigma0 p]);             % [sigma_n,p] = ExtractVectorNew(sigma0,p)     
    q = P(2:n+1)-P(1:n);
    p = p-q;                            % extracted vector
    tau = P(n+1)-sigma0;                % tau = sigma_n-sigma0 exact
    t = tp;
    tp = t + tau;                       % s = t + tau + sum(p)
    if tp==0                            % check for zero t+tau
      res = FastAccSum(p(p~=0));        % recursive call, zeros eliminated
      if rndold, setround(rndold), end
      return
    end
    q = sigma0/c3;
    u = abs(q/c4 - q);                  % u = ufp(sigma0)
    Phi = ( c5*u ) / c6;
    T = min( c7*sigma0 , c8*u );        % sum(abs(p)) <= T
    if ( abs(tp)>=Phi ) | ( 4*T<= c9 )
      tau2 = (t-tp) + tau;              % [tp,tau2] = FastTwoSum(t,tau)
      res = tp + ( tau2 + sum(p) );     % faithful.y rounded result
      if rndold, setround(rndold), end
      return
    end
  end
  