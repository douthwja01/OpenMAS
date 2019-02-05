function r = polyval(p,x,varargin)
%POLYVAL      Evaluation of polynomial p at x, same as p{x}
%
%   r = polyval(p,x)   or   r = polyval(p,x1,...,xn)
%
%A univariate polynomial p is evaluated at x.
%  For non-interval x, interval routine polyval is used.
%  For interval input x, evaluation uses either Horner's scheme or
%  sum(p(i).*x^i), the latter being faster but less accurate.
%  To switch evaluation scheme, see polynominit.
%For vector/matrix input x, result is the vector/matrix of polynomial values.
%
%For multivariate p, input x is a vector of values of the variables. Here
%  x(i) is the value for p.v{i}, i.e. evaluation depends on the order of
%  variables. To change this order, use permvars(p,perm).
%For matrix input x, rows of x are interpreted as values of the variables and
%  the result is the column vector of polynomial values.
%For multivariate polynomial p, input x may also be a cell array of values for
%  the independent variables. In that case those unknowns remain for which the 
%  empty set [] is provided and the result is a polynomial in fewer variables.
%  For example, for  p = 3 x^2 y + 4 x  , polyval(p,{5 []}) yields 75 y + 4.
%All polyval(p,1,-3), polyval(p,{1 -3}), p{1,-3} and p{{1 -3}} are the same.
%

% written  09/14/00     S.M. Rump
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

  if size(p.e,2)==1                         % univariate case
    if isa(p.c,'intval') | isa(x,'intval')
      INTLAB_POLYNOM_UEVAL = getappdata(0,'INTLAB_POLYNOM_UEVAL');
      r = intval(p.c(1))*ones(size(x));     % make sure constant polynomial produces y of same size as x
      if INTLAB_POLYNOM_UEVAL               % use Horner's scheme
        for i=2:length(p.c)
          r = r .* x + p.c(i);
        end
      else                                  % use sums of scaled powers
        n = p.e;
        if n~=0
          sizex = prod(size(x));
          if isreal(x)                        % input interval x real
            % build xpow with xpow(i) = x^i, 1<=i<=n
            if isa(x,'intval')
              xinf = x.inf(:);
              xsup = x.sup(:);
            else
              xinf = x(:);
              xsup = x(:);
            end
            index0 = ( xinf.*xsup<0 );
            xpowinf = ones(sizex,n);
            xpowsup = ones(sizex,n);
            indexposneg = ~index0;
            if any(indexposneg)                % intervals not containing zero
              k = 1;
              ypowinf = ones(sum(indexposneg),n);
              ypowsup = ypowinf;
              ypowinf(:,1) = xinf(indexposneg);
              ypowsup(:,1) = xsup(indexposneg);
              indexneg = ( ypowsup(:,1) < 0 ); % index w.r.t. ypow, not xpow
              if any(indexneg)
                dummy = ypowinf(indexneg,1);
                ypowinf(indexneg,1) = -ypowsup(indexneg,1);
                ypowsup(indexneg,1) = -dummy;
              end
              while 2*k<=n
                setround(-1)
                ypowinf(:,k+1:2*k) = (ypowinf(:,k)*ones(1,k)).*ypowinf(:,1:k);
                setround(1)
                ypowsup(:,k+1:2*k) = (ypowsup(:,k)*ones(1,k)).*ypowsup(:,1:k);
                k = 2*k;
              end
              setround(-1)
              ypowinf(:,k+1:n) = (ypowinf(:,k)*ones(1,n-k)).*ypowinf(:,1:n-k);
              setround(1)
              ypowsup(:,k+1:n) = (ypowsup(:,k)*ones(1,n-k)).*ypowsup(:,1:n-k);
              if any(indexneg)
                index = 1:2:n;
                dummy = ypowinf(indexneg,index);
                ypowinf(indexneg,index) = -ypowsup(indexneg,index);
                ypowsup(indexneg,index) = -dummy;
              end
              xpowinf(indexposneg,:) = ypowinf;
              xpowsup(indexposneg,:) = ypowsup;
            end
            if any(index0)                     % zero intervals
              k = 1;
              ypowinf = ones(sum(index0),n);
              ypowsup = ypowinf;
              ypowinf(:,1) = -xinf(index0);
              ypowsup(:,1) = xsup(index0);
              setround(1);
              while 2*k<=n
                ypowinf(:,k+1:2*k) = (ypowinf(:,k)*ones(1,k)).*ypowinf(:,1:k);
                ypowsup(:,k+1:2*k) = (ypowsup(:,k)*ones(1,k)).*ypowsup(:,1:k);
                k = 2*k;
              end
              ypowinf(:,k+1:n) = (ypowinf(:,k)*ones(1,n-k)).*ypowinf(:,1:n-k);
              ypowsup(:,k+1:n) = (ypowsup(:,k)*ones(1,n-k)).*ypowsup(:,1:n-k);
              index = ( 2:2:n );
              ypowsup(:,index) = max( ypowinf(:,index) , ypowsup(:,index) );
              ypowinf = -ypowinf;
              ypowinf(:,index) = 0;
              xpowinf(index0,:) = ypowinf;
              xpowsup(index0,:) = ypowsup;
            end
            xpow = infsup(xpowinf,xpowsup);
          else                                 % input interval x complex
            % build xpow with xpow(i) = x^i, 1<=i<=n
            if n==0                            % constant polynomial
              r = repmat(p.c,size(x));  
              setround(rndold)
              return
            else
              xpow = cintval(ones(sizex,n));
              k = 1;
              xpow(:,1) = x(:);
              while 2*k<=n
                xpow(:,k+1:2*k) = (xpow(:,k)*ones(1,k)).*xpow(:,1:k);
                k = 2*k;
              end
              xpow(:,k+1:n) = (xpow(:,k)*ones(1,n-k)).*xpow(:,1:n-k);
            end
          end
          r = sum( (ones(sizex,1)*fliplr(p.c(1:n))) .* xpow , 2 ) + p.c(n+1);
          r = reshape(r,size(x));
        end
      end
    else
      r = polyval(p.c,x);
    end
  else                                      % multivariate case
    if iscell(x)                            % x is cell array (for subpolynomials by empty components)
      k = length(x);                        % number of variables
      if k~=size(p.e,2)
        error('size of input does not match number of variables')
      end
      indexempty = logical(zeros(1,k));
      X = zeros(1,k);
      for i=1:k
        if isempty(x{i})
          indexempty(i) = 1;
        else
          if isa(x{i},'intval')
            X = intval(X);
          end
          X(i) = x{i};
        end  
      end
      if ~any(indexempty)                     % no empty components
        r = polyval(p,X);  
        setround(rndold)
        return
      end
      if all(indexempty)                      % entire polynomial
        r = p;  
        setround(rndold)
        return
      end
      r = p;
      F = (ones(size(p.e,1),1)*X(~indexempty)).^p.e(:,~indexempty);
      F(p.e(:,~indexempty)==0) = 1;
      r.c = r.c .* prod( F , 2 );
      r.e = r.e(:,indexempty);
      r.v = r.v(indexempty);
      [r.e,r.c] = collect(r.e,r.c);
      if size(r.e,2)==1                       % result univariate
        [r.e,r.c,r.v] = tounivariate(r.e,r.c,r.v);
      end
    else                                      % x is vector or matrix input
      if nargin>2                             % polyval(p,x1,x2,...)
        for i=1:length(varargin)
          if size(x,1)~=size(varargin{1},1)
            error('improper input arguments')
          end
          x = [x varargin{i}];
        end
      end
      if size(x,2)==1                         % input column vector
        x = x';                               % force row vector input
      end
      n = size(x,2);
      if n~=length(p.v)
        error('length of x does not correspond to number of variables')
      end    
      INTLAB_POLYNOM_MEVAL = getappdata(0,'INTLAB_POLYNOM_MEVAL');
      if INTLAB_POLYNOM_MEVAL                   % use Horner's scheme
        x1 = x(:,1);                            % first variable
        if n>2                                  % coefficients multivariate polynomials   
          xx = x(:,2:n);                        % remaining variables
          coeff = p;
          coeff.v = coeff.v(2:n);
          r = 0;
          for i=max(p.e(:,1)):-1:0
            index = find( p.e(:,1)==i );        % index vector to x1^i
            if any(index)
              coeff.e = p.e(index,2:n);
              coeff.c = p.c(index);           
              r = r + polyval(coeff,xx);        % update with coefficient to x1^i
            end
            if i~=0
              r = r .* x1;
            end
          end
        else                                    % coefficients univariate polynomials
          xx = x(:,2);                          % remaining variable
          coeff = p;
          coeff.v = coeff.v{2};
          r = 0;
          for i=max(p.e(:,1)):-1:0
            index = find( p.e(:,1)==i );        % index vector to x1^i
            if any(index)         
              m = max(p.e(index,2));            % degree of coefficient polynomial
              coeff.e = m;
              coeff.c = zeros(1,m+1);       
              if isa(p.c,'intval')
                coeff.c = intval(coeff.c);
              end
              coeff.c(m-p.e(index,2)+1) = p.c(index);    
              r = r + polyval(coeff,xx);        % update with coefficient to x1^i
            end
            if i~=0
              r = r .* x1;
            end                               
          end
        end
      else                                      % use sums of scaled powers
        r = 0;     
        k = size(x,1);
        if k==1                               % input one row vector
          x = ones(size(p.e,1),1)*x;
          F = x.^p.e;
          F(p.e==0) = 1;
          r = sum( p.c .* prod( F , 2) , 1 );
        else                                  % matrix input        
          maxe = max(p.e,[],1);               % maximum exponents per unknown
          powers = typeadj(zeros(k,max(maxe)+1),typeof(x)); 
          powers(:,1) = 1;
          t = ones(k,1) * ( p.c.' );
          for i=1:size(p.e,2)
            powers(:,2) = x(:,i);
            for j=2:maxe(i)
              j2 = floor(j/2);
              if 2*j2==j
                powers(:,j+1) = sqr(powers(:,j2+1));
              else
                powers(:,j+1) = powers(:,j2+1).*powers(:,j2+2);
              end
            end
            t = t .* powers(:,p.e(:,i)+1);
          end
          r = sum(t,2);
        end
      end
    end
  end
  
  setround(rndold)
