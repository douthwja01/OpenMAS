function c = vertcat(varargin)
%VERTCAT      Implements  [a(1) ; a(2) ; ...]  for intervals
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/06/04     S.M. Rump  sparse input
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
%

  a = intval(varargin{1});
  c = a;
  if a.complex
    cthin = isequal(c.rad,0);
    listcomplex = 1;
  else
    listcomplex = 0;
  end

  for i=2:length(varargin)
    a = intval(varargin{i});
    if a.complex
      if ~listcomplex
        c.complex = 1;
        e = 1e-30;
        if 1+e==1-e                           % fast check for rounding to nearest
          rndold = 0;
        else
          rndold = getround;
        end
        setround(1)
        c.rad = 0.5*(c.sup-c.inf);
        c.mid = c.inf + c.rad;
        setround(rndold)
        c.inf = [];
        c.sup = [];
        listcomplex = 1;
        cthin = 0;
      end
    end
    if listcomplex
      if ~a.complex
        a = cintval(a);
      end
      athin = isequal(a.rad,0);
      if cthin
        if ~athin
          if issparse(c.mid)
            sizecmid = size(c.mid);
            c.rad = [ sparse([],[],[],sizecmid(1),sizecmid(2),0) ; a.rad ];
          else
            c.rad = [ zeros(size(c.mid)) ; a.rad ];   % use current c.mid
          end
          cthin = 0;
        end
      else
        if athin
          if issparse(a.mid)
            sizeamid = size(a.mid);
            c.rad = [ c.rad ; sparse([],[],[],sizeamid(1),sizeamid(2),0) ];
          else
            c.rad = [ c.rad ; zeros(size(a.mid)) ];
          end
        else
          c.rad = [ c.rad ; a.rad ];
        end
      end
      c.mid = [ c.mid ; a.mid ];                    % new c.mid
    else
      c.inf = [ c.inf ; a.inf ];
      c.sup = [ c.sup ; a.sup ];
    end
  end
