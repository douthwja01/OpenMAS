function u = plus(a,b)
%PLUS         slope addition  a + b
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
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

  na = size(a.r.inf,1);
  nb = size(b.r.inf,1);

  if ( na==1 ) & ( nb~=1 )
    a.r = a.r(ones(nb,1),:);
    a.s = a.s(ones(nb,1),:);
    u.size = b.size;
  elseif ( na~=1 ) & ( nb==1 )
    b.r = b.r(ones(na,1),:);
    b.s = b.s(ones(na,1),:);
    u.size = a.size;
  else
    if ~isequal(a.size,b.size)
      error('dimensions not compatible for minus')
    end
    u.size = a.size;
  end

  u.r = a.r + b.r;
  u.s = a.s + b.s;

  u.r = rangeimprove(u);

  u = class(u,'slope');
  
  setround(rndold)

