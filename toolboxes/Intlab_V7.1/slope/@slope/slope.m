function u = slope(a,str)
%SLOPE        Slope class constructor
%
%  u = slope(a)
%
%An explicit call of the constructor is only necessary to initialize
%  a constant to be of type slope. Otherwise, any operation
%  with a dependent variable produces a result of type slope.
%
%For more details try
%
%  help slopeinit
%

%Internal representation:
%For k=1:n+1
%  u(X_1,...,X_k-1,xs_k,...,xs_n) in u.r(:,k) ,
%and for all xp_i in X_i, 1<=i<=n, there exists s in u.s with
%    u(xp) = u(xs) + u.s*(xp-xs) .
%u.r and u.s always stored as column vector; true size in u.size
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  superiorto('intval');

  if nargin==0
    u.size = [0 0];
    u.r = [];
    u.s = [];
    u = class(u,'slope');
    return
  end

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  if INTLAB_SLOPE.NUMVAR==0
    error('no slope expansion initialized')
  end

  if nargin==1

    if isa(a,'slope')
      u = a;                               % input already slope
    else
      if isa(a,'double')                   % call by slope: input constant
        a = intval(a);
      end
      if ~isa(a,'intval')
        error('invalid input for slope: input must be double or intval')
      end
      u.size = size(a);
      a = intval(a(:));
      u.r = a(:,ones(1,1+INTLAB_SLOPE.NUMVAR));
      if issparse(a)
        u.s = intval(sparse([],[],[],size(u.r,1),INTLAB_SLOPE.NUMVAR,0));
      else
        u.s = intval(zeros(size(u.r,1),INTLAB_SLOPE.NUMVAR));
      end
      u = class(u,'slope');
    end

  elseif nargin==2

    if isequal(str,'slopeinit')            % call by slopeinit

      u.size = size(a.xs);
      a.xs = a.xs(:);
      a.X = a.X(:);
      xsmat = a.xs(:,ones(1,1+INTLAB_SLOPE.NUMVAR));
      Xmat = a.X(:,ones(1,1+INTLAB_SLOPE.NUMVAR));
      u.r = tril(xsmat) + triu(Xmat,1);
      if issparse(a.xs)
        u.r = sparse(u.r);
        u.s = intval(speye(INTLAB_SLOPE.NUMVAR));
      else
        u.s = intval(eye(INTLAB_SLOPE.NUMVAR));
      end
      u = class(u,'slope');

    elseif isequal(str,'slope')            % call by @intval\slope

      a = a.init;
      u.size = size(a);
      a = intval(a(:));
      u.r = a(:,ones(1,1+INTLAB_SLOPE.NUMVAR));
      if issparse(a)
        u.s = intval(sparse([],[],[],size(u.r,1),INTLAB_SLOPE.NUMVAR,0));
      else
        u.s = intval(zeros(size(u.r,1),INTLAB_SLOPE.NUMVAR));
      end
      u = class(u,'slope');

    else
      error('invalid call of slope constructor')
    end

  else

    error('invalid call of slope constructor')
  end

  % avoid Matlab 6.5f bug: 
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if prod(size(u.r))==1
    u.r = full(u.r);
  end
  