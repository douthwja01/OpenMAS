function r = hessian(a,str)
%HESSIAN      Hessian class constructor
%
%  r = hessian(a)
%
%An explicit call of the constructor is only necessary to initialize
%  a constant to be of type hessian. Otherwise, any operation
%  with a dependent variable produces a result of type hessian.
%
%For more details see
%
%  help hessianinit
%

%
%For N (=INTLAB_HESSIAN_NUMVAR) denoting the number of independent variables, storage conventions are as follows:
%  r.x    is scalar, vector or matrix; size(r) = size(r.x)
%  r.dx   is the gradient of transpose(r(:)), that is the matrix so that column i, 
%           i.e. r.dx(:,i), of length N being equal to the gradient of r(i).
%           ===> Careful, columns of r.dx are gradients <===
%  r.hx   stores information for the Hessian of transpose(r(:)). Storage conventions are as follows:
%           Let H := r.hx(:,i), then reshape(H+H.',N,N) is the true Hessian matrix of r(i).
%           ===> Careful, columns of r.hx are a sort of unsymmetric part of Hessian <===
%         For sparse matrices this storage scheme may reduce memory requirements significantly.
%

% written  04/04/04     S.M. Rump  Thanks to Arnold for suggesting .hx, is nicer than previous .ddx
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  superiorto('intval');

  if nargin==0
    r.x = [];
    r.dx = [];
    r.hx = [];
    r = class(r,'hessian');
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  if N==0
    error('no dependent variables initialized for use of hessian')
  end

  if nargin==1

    if isa(a,'hessian')
      r = a;
    else
      r.x = a;
      len = prod(size(r.x));
      INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
      if N<INTLAB_HESSIAN_SPARSE
        r.dx = zeros(N,len);
        r.hx = zeros(N*N,len);
      else
        r.dx = sparse([],[],[],N,len);
        r.hx = sparse([],[],[],N*N,len);
      end
      r = class(r,'hessian');
    end

  elseif nargin==2

    if isequal(str,'hessianinit')         % call by hessianinit

      r.x = a.init;
      INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
      if N<INTLAB_HESSIAN_SPARSE
        r.dx = eye(N);
        r.hx = zeros(N*N,N);
      else
        r.dx = speye(N);
        r.hx = sparse([],[],[],N*N,N);
      end
      if isa(r.x,'intval')
        r.dx = intval(r.dx);
        r.hx = intval(r.hx);
      end
      r = class(r,'hessian');

    elseif isequal(str,'hessian')         % call by @intval\hessian

      r.x = a.init;
      len = prod(size(r.x));
      INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
      if N<INTLAB_HESSIAN_SPARSE
        r.dx = intval(zeros(N,len));
        r.hx = intval(zeros(N*N,len));
      else
        r.dx = intval(sparse([],[],[],N,len));
        r.hx = intval(sparse([],[],[],N*N,len));
      end
      r = class(r,'hessian');

    elseif isequal(str,'hessianintval')   % call by @intval\intval

      r.x = intval(a.x);
      r.dx = intval(a.dx);
      r.hx = intval(a.hx);
      r = class(r,'hessian');

    elseif isequal(str,'random')          % generates .dx and .hx randomly, only for test purposes

      if isa(a,'struct')                  % input interval
        a = a.init;
      end
      r.x = a;
      len = prod(size(r.x));
      INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
      if N<INTLAB_HESSIAN_SPARSE
        if isreal(r.x)
          r.dx = random(N,len);
          r.hx = random(N^2,len);
        else
          r.dx = randomc(N,len);
          r.hx = randomc(N^2,len);
        end
      else
        dendx = max(5/N,1/N/len);
        denhx = max(5/(N^2),1/N^2/len);
        % about 5 nonzero derivatives per unknown
        if isreal(r.x)
          r.dx = 2*sprand(N,len,dendx);
          r.dx = r.dx - spones(r.dx);
          r.hx = 2*sprand(N^2,len,denhx);
          r.hx = r.hx - spones(r.hx);
        else
          r.dx = 2 * ( sprand(N,len,dendx) + sqrt(-1)*sprand(N,len,dendx) );
          r.dx = r.dx - (1+sqrt(-1))*spones(r.dx);
          r.hx = 2 * ( sprand(N^2,len,denhx) + sqrt(-1)*sprand(N^2,len,denhx) );
          r.hx = r.hx - (1+sqrt(-1))*spones(r.hx);
        end
      end
      if isa(r.x,'intval')
        r.dx = intval(r.dx);
        r.hx = intval(r.hx);
      end
      r = class(r,'hessian');

    else
      error('invalid call of hessian constructor')
    end

  else

    error('invalid call of constructor hessian')
  end
  
  % avoid Matlab 6.5f bug: 
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if prod(size(r.x))==1
    r.x = full(r.x);
  end
  if prod(size(r.dx))==1
    r.dx = full(r.dx);
  end
  if prod(size(r.hx))==1
    r.hx = full(r.hx);
  end
    
  if rndold
    setround(rndold)
  end
