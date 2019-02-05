function c = minus(a,b)
%MINUS        Implements  a - b  for intervals
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    remove check for 'double'
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(a,'intval')                   % a is double
    if ~isreal(a) | b.complex           % complex case
      if ~b.complex                     % C - IR
        setround(1)
        b.mid = b.inf + 0.5*(b.sup-b.inf);
        b.rad = b.mid - b.inf;
      end
      c.complex = 1;                    % R - IC  or  C - IC
      c.inf = [];
      c.sup = [];
      setround(-1)
      c1 = a - b.mid;
      setround(1)
      c.mid = a - b.mid;
      if isequal(b.rad,0)
        c.rad = abs(c.mid-c1);
      else
        c.rad = abs(c.mid-c1) + b.rad;
      end
      if isequal(c.rad,0)
        c.rad = 0;
      end
    else                                  % real case  R - IR
      c.complex = 0;
      setround(-1)
      c.inf = a - b.sup;
      setround(1)
      c.sup = a - b.inf;
      c.mid = [];
      c.rad = [];
    end
  elseif ~isa(b,'intval')                 % b is double
    if a.complex | ~isreal(b)             % complex case
      if ~a.complex                       % IR - C
        setround(1)
        a.mid = a.inf + 0.5*(a.sup-a.inf);
        a.rad = a.mid - a.inf;
      end
      c.complex = 1;                      % IC - R  or  IC - C
      c.inf = [];
      c.sup = [];
      setround(-1)
      c1 = a.mid - b;
      setround(1)
      c.mid = a.mid - b;
      if isequal(a.rad,0)
        c.rad = abs(c.mid-c1);
      else
        c.rad = abs(c.mid-c1) + a.rad;
      end
      if isequal(c.rad,0)
        c.rad = 0;
      end
    else                                  % real case  IR - R
      c.complex = 0;
      setround(-1)
      c.inf = a.inf - b;
      setround(1)
      c.sup = a.sup - b;
      c.mid = [];
      c.rad = [];
    end
  else                                    % both a and b interval
    if a.complex | b.complex              % complex case
      if ~a.complex
        setround(1)
        a.mid = a.inf + 0.5*(a.sup-a.inf);
        a.rad = a.mid - a.inf;
      end
      if ~b.complex
        setround(1)
        b.mid = b.inf + 0.5*(b.sup-b.inf);
        b.rad = b.mid - b.inf;
      end
      c.complex = 1;                      % IC - IC
      c.inf = [];
      c.sup = [];
      setround(-1)
      c1 = a.mid - b.mid;
      setround(1)
      c.mid = a.mid - b.mid;
      c.rad = abs(c.mid-c1);
      if ~isequal(a.rad,0)
        c.rad = c.rad + a.rad;
      end
      if ~isequal(b.rad,0)
        c.rad = c.rad + b.rad;
      end
      if isequal(c.rad,0)
        c.rad = 0;
      end
    else                                  % real case  IR - IR
      c.complex = 0;
      setround(-1)
      c.inf = a.inf - b.sup;
      setround(1)
      c.sup = a.sup - b.inf;
      c.mid = [];
      c.rad = [];
    end
  end

  c = class(c,'intval');

  setround(rndold)
