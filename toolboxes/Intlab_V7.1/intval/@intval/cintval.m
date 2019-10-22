function a = cintval(a)
%CINTVAL      Type cast to complex interval
%
%  c = cintval(a)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 09/02/00     S.M. Rump  same as tocmplx, rounding preserved
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/03/05     S.M. Rump  sparse input
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
%

  if ~a.complex
    if isequal(a.inf,a.sup)
      a.complex = 1;
      a.mid = a.inf;
      a.rad = 0;
    else
      if issparse(a.inf)                  % sparse input
        [i,j,s] = find(a);
        [m,n] = size(a.inf);
        a = sparse(i,j,cintval(full(s)),m,n);
        return
      else                                % full input
        a.complex = 1;
        index = ( a.inf==a.sup );
        anyindex = any(index(:));
        if anyindex
          a.mid = a.inf;
          a.rad = zeros(size(a.inf));
        end
        e = 1e-30;
        if 1+e==1-e                       % fast check for rounding to nearest
          rndold = 0;
        else
          rndold = getround;
        end
        if ~anyindex
          setround(1)
          a.mid = a.inf + 0.5*(a.sup-a.inf);
          a.rad = a.mid - a.inf;
        else
          index = ~index;
          setround(1)
          ainf = a.inf(index);
          a.mid(index) = ainf + 0.5*(a.sup(index)-ainf);
          a.rad(index) = a.mid(index) - ainf;
        end
        setround(rndold)                  % set rounding to previous value
      end
    end
    a.inf = [];
    a.sup = [];
  end
