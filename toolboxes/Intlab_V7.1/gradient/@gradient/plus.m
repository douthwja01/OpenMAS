function r = plus(a,b)
%PLUS         Gradient addition  a + b
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/23/06     S.M. Rump  comments corrected
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(a,'gradient')         % non-gradient plus gradient
    r.x = a + b.x;
    if prod(size(b.x))>1        % b is not scalar
      r.dx = b.dx;
    else                        % b is scalar, a may be array
      r.dx = b.dx(ones(prod(size(r.x)),1),:);
    end
    if isa(a,'intval') & isa(b.x,'double')
      r.dx = intval(r.dx);
    end
  elseif ~isa(b,'gradient')     % gradient plus non-gradient
    r.x = a.x + b;
    if prod(size(a.x))>1        % a is not scalar
      r.dx = a.dx;
    else                        % a is scalar, b may be array
      r.dx = a.dx(ones(prod(size(r.x)),1),:);
    end
    if isa(b,'intval') & isa(a.x,'double')
      r.dx = intval(r.dx);
    end
  else                          % gradient plus gradient
    r.x = a.x + b.x;
    if prod(size(a.x))==1       % scalar gradient plus gradient
      r.dx = a.dx(ones(size(b.dx,1),1),:) + b.dx;
    elseif prod(size(b.x))==1
      r.dx = a.dx + b.dx(ones(size(a.dx,1),1),:);
    else
      r.dx = a.dx + b.dx;
    end
  end

  r = class(r,'gradient');
  
  if rndold
    setround(rndold)
  end
