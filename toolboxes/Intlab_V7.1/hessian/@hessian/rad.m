function r = rad(a)
%RAD          Hessian radius
%
%  r = rad(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if issparse(a.hx)
    if isa(a.x,'intval')
      r.x = rad(a.x);
      if isequal(r.x,0)
        [m n] = size(a.x);
        r.x = sparse([],[],[],m,n,0);
      end
      r.dx = rad(a.dx);
      if isequal(r.dx,0)
        [m n] = size(a.dx);
        r.dx = sparse([],[],[],m,n,0);
      end
      r.hx = rad(a.hx);
      if isequal(r.hx,0)
        [m n] = size(a.hx);
        r.hx = sparse([],[],[],m,n,0);
      end
    else
      [m n] = size(a.x);
      r.x = sparse([],[],[],m,n,0);
      [m n] = size(a.dx);
      r.dx = sparse([],[],[],m,n,0);
      [m n] = size(a.hx);
      r.hx = sparse([],[],[],m,n,0);
    end
  else
    if isa(a.x,'intval')
      r.x = rad(a.x);
      if isequal(r.x,0)
        r.x = zeros(size(a.x));
      end
      r.dx = rad(a.dx);
      if isequal(r.dx,0)
        r.dx = zeros(size(a.dx));
      end
      r.hx = rad(a.hx);
      if isequal(r.hx,0)
        r.hx = zeros(size(a.hx));
      end
    else
      r.x = zeros(size(a.x));
      r.dx = zeros(size(a.dx));
      r.hx = zeros(size(a.hx));
    end
  end
  
  r = class(r,'hessian');
