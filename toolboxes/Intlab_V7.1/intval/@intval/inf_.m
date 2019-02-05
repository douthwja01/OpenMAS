function c = inf_(a)
%INF_         Implements  inf(a)  (cures problems with inf)
%
%   c = inf(a)
%
% On return, inf(a) <= alpha for all entries alpha in a
%
%************************************************************************
%********  due to conflict with internal variable inf (infinity)  *******
%********                    use function inf_                    *******
%************************************************************************
%

% written  10/05/99     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/10/07     S.M. Rump  performance, huge arrays
%

  if a.complex
    if isequal(a.rad,0)                  % faster for sparse matrices
      c = a.mid;
    else
      e = 1e-30;
      if 1+e==1-e                        % fast check for rounding to nearest
        rndold = 0;
      else
        rndold = getround;
      end
      setround(-1)
      if issparse(a.rad)
        [m,n] = size(a.rad);
        [I,J,arad] = find(a.rad);
        c = a.mid - sparse(I,J,complex(arad,arad),m,n);
      else
        c = a.mid - complex(a.rad,a.rad);
      end
      setround(rndold)
    end
  else
    c = a.inf;
  end
