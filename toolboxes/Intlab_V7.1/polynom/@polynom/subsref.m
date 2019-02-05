function r = subsref(p,s)
%SUBSREF      Implements subscripted references for polynomials
%
%  r = p(i)
%
%For univariate polynomials p, p(i) is the coefficient of x^i, 0<=i<=degree(p).
%  For i>degree(p), p(i):=0. Similarly, p(i:j) is the vector of coefficients
%  [ p(i) p(i+1) ... p(j) ], or, p(:) is the (row) vector of all coefficients
%  of p, the same as vector(p) for univariate polynomials.
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
%Polynomial evaluation is denoted by p{x} or p{x1,...,xn}, computing the value 
%  of p at x. This is the same as polyval(p,x) of polyval(p,x1,...,xn).
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
%

  while 1
    if ~isa(p,'polynom')
      r = subsref(p,s(1));
    elseif strcmp(s(1).type,'()')     % index reference p(i)
      if size(p.e,2)==1               % univariate polynomial
        if length(s(1).subs)>1
          error('invalid call: more than one index')
        end
        if isequal(s(1).subs{1},':')
          r = p.c;
        else
          n1 = length(p.c);
          index = n1-s(1).subs{1};    % index vector, p(0) constant term, etc.
          if ~isreal(index) | ~isequal(index,round(index))
            error('index must be integer')
          end
          if any(index>n1)            % index negative
            error('index out of range')
          end
          if any(index<1)             % index greater than degree
            indexgt0 = ( index>0 );
            r = typeadj( zeros(1,length(index)) , typeof(p.c) );
            if any(indexgt0)
              r(indexgt0) = p.c(index(indexgt0));
            end
          else
            r = p.c(index);
          end
        end
      else                            % multivariate polynomial
        k = size(p.e,2);              % number of variables
        if length(s(1).subs)>k
          error('too many indices')
        end
        index = logical(zeros(1,k));
        exponents = [];
        for i=1:length(s(1).subs)
          if isempty(s(1).subs{i}) | isequal(s(1).subs{i},':')
            index(i) = 0;
          else
            index(i) = 1;
            exponents = [ exponents s(1).subs{i} ];
          end
        end 
        if isempty(exponents)
          r = p;
          return 
        end
        I = all( ( p.e(:,index) == repmat(exponents,size(p.e,1),1) ) , 2 );
        if ~any(I)                    % coefficients do not occur
          r = typeadj( 0 , typeof(p.c) );
          return
        end
        r.e = p.e(I,~index);
        if isempty(r.e)               % coefficient is single constant
          r = p.c(I);
          return
        end
        r.c = p.c(I);
        r.v = p.v(~index);
        if size(r.e,2)==1             % coefficient polynomial univariate
          if iscell(r.v)              % be sure r.v is char for univ. pol.
            r.v = r.v{1};
          end
          n = max(r.e);               % r.e is vector
          c = r.c;
          r.c = typeadj( zeros(1,n+1) , typeof(c) );
          r.c(n+1-r.e) = c;
          r.e = n;
          if r.e==0
            r = r.c;
            return
          end
        else
          r = normalize(r);
        end
        r = class(r,'polynom');
      end
    elseif strcmp(s(1).type,'{}')     % polynomial evaluation p{x}
      r = polyval(p,s(1).subs{:});
    elseif strcmp(s(1).type,'.')      % polynomial access to inf, sup, ...
      if     strcmp(s(1).subs,'mid'), r = mid(p);
      elseif strcmp(s(1).subs,'rad'), r = rad(p);
      elseif strcmp(s(1).subs,'inf'), r = inf(p);
      elseif strcmp(s(1).subs,'sup'), r = sup(p);
      elseif strcmp(s(1).subs,'e')
        if iscell(p.v)                  % multivariate
          r = flipud(sortrows(p.e));     
        else
          r = p.e;                      % univariate
        end
      elseif strcmp(s(1).subs,'c')
        if iscell(p.v)                  % multivariate          
          [ dummy index ] = sortrows(p.e); r = flipud(p.c(index));
        else
          r = p.c;
        end
      elseif strcmp(s(1).subs,'v'), r = p.v;
      else
        error('invalid reference for polynomial')
      end
    else
      error('invalid index reference for polynomial')
    end
    if length(s)==1
      return
    end
    s = s(2:end);
    p = r;
  end
