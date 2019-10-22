function r = rad(a)
%RAD          Gradient radius
%
%  r = rad(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isa(a.x,'intval')
    r.x = rad(a.x);
    if isequal(r.x,0)
      if issparse(a.x)
        r.x = sparse([],[],[],size(a.x,1),size(a.x,2));
      else
        r.x = zeros(size(a.x));
      end
    end
    r.dx = rad(a.dx);
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
