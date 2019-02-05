function p = subsasgn(p,s,r)
%SUBSASGN     Implements subscripted assignment for polynomials
%
%  p(i) = r
%
%For univariate polynomials p, p(i) is the coefficient of x^i, 0<=i<=degree(p).
%  For i>degree(p), p(i):=0. Similarly, p(i:j) is the vector of coefficients
%  [ p(i) p(i+1) ... p(j) ], or, p(:) is the (row) vector of all coefficients
%  of p, the same as vector(p) for univariate polynomials. Especially, 
%  p(i) = [] of p(i:j) = [] cancels the i-th or i..j-th coefficient, respectively.
%
%For multivariate polynomials p in k variables x_1..x_k, p(i) is the
%  coefficient of x_1^i, i.e. a polynomial in k-1 variables.
%  Similarly, p(i,[],j) or p(i,:,j) is the coefficient polynomial of
%  x1^i*x3^j. Indices i,j,... must be single indices, no range.
%Note that access to coefficients refers to the current order p.v of
%  variables of p. To change this order, see permvars.
%For example, for a polynomial in three variables p.v={'x','y','z'},
%  p(1,0,3) is the coefficient of x*z^3 (a constant), where p(1,[],3)
%  if the coefficient of x*z^3, a univariate polynomial in y.
%
%Polynomial evaluation is denoted by p{x}, computing the value of p at x.
%  This is the same as polyval(p,x).
%For univariate polynomials, x may be a vector or matrix yielding the vector
%  or matrix of polynomial values evaluated at the corresponding coefficients.
%For multivariate polynomials, x is a vector of values of the variables.
%  For x being a matrix, the result is the (column) vector of p{x(i,:)}.
%
%Moreover, p.mid, p.rad, p.inf, p.sup give access to the midpoint, radius,
%  infimum and supremum of p, respectively.
%
%Finally, p.e, p.c and p.v give access to the arrays of exponents, coefficients
%  and variables of p, such that polynom(p.e,p.c,p.v) is again p. Single variables
%  are accessed by p.v{i}.
%
%In the univariate case,  p == polynom(p.c,p.v)  [ degree is length(p.c)+1 ].
%In the multivariate case,  p == polynom(p.e,p.c,p.v)  [ p.e is (sparse) exponent set ].
%

% written  11/20/97     S.M. Rump
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

  while 1
    if ~isa(p,'polynom')
      p = subsasgn(p,s(1),r);
    elseif strcmp(s(1).type,'()')     % index reference p(i)
      if size(p.e,2)==1               % univariate polynomial
        if length(s(1).subs)>1
          error('invalid call: more than one index')
        end
        if isequal(s(1).subs{1},':')
          p.c(:) = typeadj(r,typeof(p.c));
        else
          n1 = length(p.c);
          index = n1-s(1).subs{1};     % index vector, p(0) constant term, etc.
          if ~isreal(index) | ~isequal(index,round(index))
            error('index must be integer')
          end
          if any(index>n1)             % index negative
            error('index out of range')
          end
          if ( length(r)>1 ) & ( length(index)~=length(r) )
            error('length of arguments do not match, invalid assignment')
          end          
          if isa(r,'intval') & ~isa(p.c,'intval')
            p.c = intval(p.c);
          end
          m = min(index);
          if isempty(r)
            r = 0;
          end
          if m<=0                      % index greater than degree
            p.c = [zeros(1,-m+1) p.c];
            p.c(index-m+1) = r;
          else
            p.c(index) = r;
          end
        end
        m = min(find(p.c~=0));
        if isempty(m)                  % zero polynomial
          p.c = typeadj(0,typeof(p.c));
        elseif m>1                     % leading coefficients zero
          p.c = p.c(m:length(p.c));    % Matlab bug: ..end does not work for user-defined data types
        end
        p.e = length(p.c)-1;
      else                             % multivariate polynomial
        k = size(p.e,2);               % number of variables
        if length(s(1).subs)>k
          error('too many indices')
        end
        exponents = zeros(1,k);
        for i=1:length(s(1).subs)
          if isempty(s(1).subs{i}) | isequal(s(1).subs{i},':')
            error('Only integer indices allowed in multivariate polynomial assignment')
          end
          exponents(i) = s(1).subs{i};
        end 
        I = all( ( p.e == repmat(exponents,size(p.e,1),1) ) , 2 );
        if any(I)                     % coefficient does occur
          if r==0
            p.e(I,:) = [];
            p.c(I) = [];
          else
            p.c(I) = r;
          end      
        else                          % coefficient does not occur
          if r~=0
            p.e = [p.e ; exponents];
            if isa(r,'intval')
              p.c = intval(p.c);
            end
            p.c = [p.c ; r];
          end
        end
        p = normalize(p);
      end
    elseif strcmp(s(1).type,'.')      % polynomial access to .v
      if strcmp(s(1).subs,'v')
        if ischar(p.v)                % univariate p
          if ~ischar(r)
            error('variable of univariate polynomial must be string')
          else
            p.v = r;
          end
        else                          % multivariate p
          if ( length(s)==2 ) & ( s(2).type=='{}' )
            if s(2).subs{1}>length(p.v)
              error('index of variable too large')
            end
            if ~ischar(r)
              error('variable name must be string')
            end
            I = find(strcmp(p.v,r));
            if ~isempty(I) & ( I~=s(2).subs{1} )
              error('duplicate variable name not allowed')
            end
            p.v{s(2).subs{1}} = r;  
            if rndold
              setround(rndold)
            end
            return
          elseif ~iscell(r)
            error('variables of multivariate polynomials must be cell array of strings')
          else
            if length(r)~=length(p.v)
              error('number of variables does not match')
            else
              p.v = r;
            end
          end
        end
      else
        error('invalid reference for polynomial')
      end
    else
      error('invalid index reference for polynomial')
    end
    if length(s)==1  
      if rndold
        setround(rndold)
      end
      return
    end
    error('invalid call of polynom/subsasgn')
  end
