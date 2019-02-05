function r = polynom(pc,pe,pv)
%POLYNOM      Polynomial class constructor
%
%Polynomial coefficients maybe real, complex, or intval
%There are various ways to generate a polynomial:
%
%  p = polynom(pc)         Univariate polynomial to variable 'x', where pc is
%                            array of n+1 coefficients p_n, ... , p_1 , p_0
%
%  p = polynom(pc,var)     As above but to variable var, where var is a string
%
%  p = polynom(pc,pe)      Univariate polynomial to variable 'x' with coefficients pc
%                            to exponents pe
%
%  p = polynom(pc,pe,var)  As above but to variable var, where var is a string
%
%  p = polynom(pc,pe,pv)   If pv is cell array of strings, then p is multivariate polynomial with
%                            pc  m x 1 array of corresponding coefficients
%                            pe  m x n array of nonnegative integer exponents
%                            pv  cell array with n components for independent variables
%                          Result p is univariate in case pv consists of only one variable.
%
%There is access to coefficients, exponents and variables:
%For univariate p          p.c is the 1 x (n+1) array of coefficients
%                          p.e is the degree n (0 for constant polynomial)
%                          p.v is the string of the independent variable
%
%For multivariate p        p.c is the m x 1 array of coefficients
%                          p.e is the m x n array of exponents
%                          p.v is the n cell array of strings of the independent variables
%
%Examples:
%
%  r = polynom([1 0 -2])
%  r = polynom([1 -2],[2 0])
%  r = polynom([1 0 -2],'x')
%  r = polynom([1 -2],[2 0],'x')
%  pc = [ 1 ; -2 ],  pe = [ 2 ; 0 ],  pv = {'x'},  r = polynom(pc,pe,pv)
%
%all generate the same polynomial  x^2-2 .
%
%  pc = [ 1 ; 2 ; 1 ],  pe = [ 2 0 ; 1 1 ; 0 2 ],  pv = {'x','y'},
%  r = polynom(pc,pe,pv)
%
%generates  x^2 + 2xy + y^2 .
%
%For access to coefficients like r(0), evaluation r{x} at a point x,
%  cf.  help polynom\subsref
%

%p is univariate iff ischar(p.v) (iscell for multivariate polynomial),
%  and iff size(p.e,2)==1 (=n # of variables for multivariate polynomial)
%  Variable string may be empty (by removevars) for constant polynomial. In this
%  polynomial is compatible with every other polynomial.
%
%Class polynom:  univariate case
%  p.e   degree n of polynomial
%  p.c   array (row) of n+1 coefficients p_n, ..., p_1, p_0
%  p.v   string of variable name (default 'x')
%
%Class polynom:  multivariate case
%  p.e   m x n integer array of exponents, m coefficients to n variables
%  p.c   array of m coefficients corresponding to p.e
%          suitable types: double, complex, intval
%          zero components are eliminated
%  p.v   cell array of n variable names (default 'x1', 'x2', ...)
%

% written  09/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/11/07     S.M. Rump  typo
%

  superiorto('intval');
  
  if nargin==0                          % call r = polynom
    r.e = [];
    r.c = [];
    r.v = [];
    r = class(r,'polynom');
  elseif nargin==1                      % call r = polynom(pc)
    if isa(pc,'polynom')
      r = pc;
    else
      [m n] = size(pc);
      if ( m~=1 ) & ( n~=1 )
        error('invalid call of constructor polynom: input must be vector')
      end
      r.e = length(pc)-1;
      r.c = pc(:).';
      r.v = 'x';
      r = normalize( class(r,'polynom') );
    end
  elseif nargin==2
    if ischar(pe)                       % call r = polynom(pc,var)   [var=pe]
      [m n] = size(pc);
      if ( m~=1 ) & ( n~=1 )
        error('invalid call of constructor polynom: input must be vector')
      end
      r.e = length(pc)-1;
      r.c = pc(:).';    
      r.v = pe;
      if ( r.e~=0 ) & isempty(r.v)
        error('no variable only allowed for constants')
      end
    else                                % call r = polynom(pc,pe)
      if ~isequal(size(pc),size(pe))
        error('sizes of coefficients and exponents do not match')
      end
      if ( size(pc,1)~=1 ) & ( size(pc,2)~=1 )
        error('arrays of exponents and coefficients of univariate polynomials must be vectors')
      end
      if ~isequal(pe,round(pe))
        error('degree must be integer')
      end
      r.e = max(pe);
      r.c = typeadj(zeros(1,r.e+1),typeof(pc));
      r.c(r.e-pe+1) = pc;
      r.v = 'x';
    end
    r = normalize( class(r,'polynom') );
  elseif nargin==3                      % call r = polynom(pc,pe,pv)
    if iscell(pv) & ( length(pv)==1 )   % univariate polynomial
      pv = pv{1};
    end
    if ischar(pv)                       % univariate polynomial r = polynom(pc,pe,var)   [var=pv]
      if ~isequal(size(pc),size(pe))
        error('sizes of exponents and coefficients do not match')
      end
      if ( size(pc,1)~=1 ) & ( size(pc,2)~=1 )
        error('arrays of exponents and coefficients of univariate polynomials must be vectors')
      end
      if ~isequal(pe,round(pe))
        error('degree must be integer')
      end
      r.e = max(pe);
      r.c = typeadj(zeros(1,r.e+1),typeof(pc));
      r.c(r.e-pe+1) = pc;
      r.v = pv;
      if ( r.e~=0 ) & isempty(r.v)
        error('polynomial without dependent variable only allowed for constants')
      end
    else                                % multivariate polynomial polynom(pc,pe,cell)
      r.e = pe;
      r.c = pc;
      r.v = pv;
      [m n] = size(r.e);
      if size(r.c,1)==1
        r.c = r.c.';
      end
      if ~isequal(size(r.c),[m 1])
        error('sizes of arrays for exponents and coefficients do not match')
      end
      if length(r.v)~=n
        error('sizes of arrays for exponents and variables do not match')
      end
      for i=1:n
        if ~ischar(r.v{i}) | isempty(r.v{i})
          error('independent variables must be strings')
        end
      end
    end
    r = normalize( class(r,'polynom') );
  else
    error('invalid call of constructor polynom')
  end

  % avoid Matlab 6.5f bug: 
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if prod(size(r.c))==1
    r.c = full(r.c);
  end
