function r = diam(a)
%DIAM         Gradient diameter
%
%  r = diam(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
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

  if isa(a.x,'intval')
    r.x = diam(a.x);
    if isequal(r.x,0)
      if issparse(a.x)
        r.x = sparse([],[],[],size(a.x,1),size(a.x,2));
      else
        r.x = zeros(size(a.x));
      end
    end
    r.dx = diam(a.dx);
    if isequal(r.dx,0)
      if issparse(a.dx)
        r.dx = sparse([],[],[],size(a.dx,1),size(a.dx,2));
      else
        r.dx = zeros(size(a.dx));
      end
    end
    r = class(r,'gradient');
  else
    if issparse(a.x)
      r.x = sparse([],[],[],size(a.x,1),size(a.x,2));
    else
      r.x = zeros(size(a.x));
    end
    if issparse(a.dx)
      r.dx = sparse([],[],[],size(a.dx,1),size(a.dx,2));
    else
      r.dx = zeros(size(a.dx));
    end
    r = class(r,'gradient');
  end
  
  if rndold
    setround(rndold)
  end
