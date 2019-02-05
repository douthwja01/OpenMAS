function u = mtimes(a,b)
%MTIMES       Slope multiplication  a * b
%
%One argument must be scalar.
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input, improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  if ~isa(a,'slope')
    a = slope(a);
  end
  if ~isa(b,'slope')
    b = slope(b);
  end

  na = size(a.r,1);
  nb = size(b.r,1);
  scalar = 0;

  if ( na==1 ) & ( nb~=1 )
    a.r = a.r(ones(nb,1),:);
    a.s = a.s(ones(nb,1),:);
    na = nb;
    u.size = b.size;
    scalar = 1;
  elseif ( na~=1 ) & ( nb==1 )
    b.r = b.r(ones(na,1),:);
    b.s = b.s(ones(na,1),:);
    nb = na;
    u.size = a.size;
    scalar = 1;
  elseif ( length(a.size)>2 ) | ( length(b.size)>2 ) | ...
           ~isequal(a.size(2),b.size(1))
    error('Dimension of factors in slope/mtimes not compatible')
  end

  if scalar                 % one factor scalar
    u.r = a.r .* b.r;
    indexc = 1:INTLAB_SLOPE.NUMVAR;
    indexr = 2:INTLAB_SLOPE.NUMVAR+1;
    u.s = intersect( a.r(:,indexr) .* b.s + b.r(:,indexc) .* a.s , ...
                     b.r(:,indexr) .* a.s + a.r(:,indexc) .* b.s );
  else                      % both factors non-scalar
    u.size = [ a.size(1) b.size(2) ];
    if issparse(a.s) & issparse(b.s)
      u.r = intval(sparse([],[],[],prod(u.size),INTLAB_SLOPE.NUMVAR+1,0));
      u.s = intval(sparse([],[],[],prod(u.size),INTLAB_SLOPE.NUMVAR,0));
    else
      u.r = intval(zeros([prod(u.size) INTLAB_SLOPE.NUMVAR+1]));
      u.s = intval(zeros([prod(u.size) INTLAB_SLOPE.NUMVAR]));
    end
    for i=1:INTLAB_SLOPE.NUMVAR+1
      uri = reshape(a.r(:,i),a.size) * reshape(b.r(:,i),b.size);
      u.r(:,i) = uri(:);
    end
    for i=1:INTLAB_SLOPE.NUMVAR
      f1 = reshape(a.r(:,i+1),a.size) * reshape(b.s(:,i),b.size) + ...
           reshape(a.s(:,i),a.size) * reshape(b.r(:,i),b.size) ;
      % use  a*b = (b'*a')'
      f2 = ( reshape(b.r(:,i+1),b.size)' * reshape(a.s(:,i),a.size)' + ...
             reshape(b.s(:,i),b.size)' * reshape(a.r(:,i),a.size)' )' ;
      usi = intersect(f1(:),f2(:));
      u.s(:,i) = usi(:);
    end
  end

  u.r = rangeimprove(u);

  u = class(u,'slope');
  
  setround(rndold)
