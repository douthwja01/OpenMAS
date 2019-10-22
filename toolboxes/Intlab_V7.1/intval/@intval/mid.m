function c = mid(a)
%MID          Implements  mid(a)  for intervals (rounded)
%
%   c = mid(a)
%
% mid(a) and rad(a) computed such that
%    alpha  in  < mid(a) , rad(a) >  for all alpha in a
%
%For intervals at least 3 bits wide, the midpoint is always an inner point.
%

% written  10/16/98     S.M. Rump
% modified 06/22/99     S.M. Rump  for sparse matrices
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 09/10/07     S.M. Rump  performance, huge arrays
% modified 10/18/08     S.M. Rump  again huge arrays
% modified 07/23/09     S.M. Rump  changed formula: midpoint now in rnd to nearest
%                                    to make sure mid([1-eps,1+2eps])=1 
%                                    (thanks to Gerhard Heindl for pointing to this)
% modified 10/06/09     S.M. Rump  check for rndold
%

  if a.complex                            % complex or thin interval
    c = a.mid;
  else
    e = 1e-30;
    if 1+e==1-e                         % fast check for rounding to nearest
      rndold = 0;
    else
      rndold = getround;
      setround(0)
    end
    % use a.inf + (0.5*a.sup-0.5*a.inf) for correct result in case a.sup-a.inf overflows
    [m,n] = size(a.inf);
    if m*n<2^31                           % input not huge
      c = a.inf + (0.5*a.sup-0.5*a.inf);  % make sure result correct in underflow range
      indexinf = isinf(a.inf);
      indexsup = isinf(a.sup);
      anyindexinf = any(indexinf(:));
      anyindexsup = any(indexsup(:));
      if anyindexinf                      % make sure mid([-inf,x]) is x
        c(indexinf) = a.sup(indexinf);
      end
      if anyindexsup                      % make sure mid([x,inf]) is x
          c(indexsup) = a.inf(indexsup);
      end
      if anyindexinf | anyindexsup        % some components are inf
        c(indexinf & indexsup) = 0;       % make sure mid([-inf,inf]) is 0
      end
    else                                  % take care of huge matrices
      % check some components are inf
      % careful with intervals [0,2] or [-2,0] or [-2,2]
      [Iinf,Jinf,Sinf] = find(a.inf);
      [Isup,Jsup,Ssup] = find(a.sup);
      if ( ~isempty(Iinf) ) | ( ~isempty(Isup) )
        ainfsup = sparse([Iinf;Isup],[Jinf;Jsup],[complex(Sinf,0);complex(0,Ssup)],m,n);
        [I,J,S] = find(ainfsup);
        Sold = S;
        S = real(S) + (0.5*imag(S)-0.5*real(S));
        indexinf = isinf(real(Sold));
        if any(indexinf(:))               % make sure mid([-inf,x]) is x
          S(indexinf) = imag(Sold(indexinf));
        end
        indexsup = isinf(imag(Sold));
        if any(indexsup(:))               % make sure mid([x,inf]) is x
          S(indexsup) = real(Sold(indexsup));
        end
        S(indexinf & indexsup) = 0;       % make sure mid([-inf,inf]) is 0
        c = sparse(I,J,S,m,n);
      else
        c = a.inf + 0.5*(a.sup-a.inf);    % make sure result correct in underflow range
      end
    end
    if rndold
      setround(rndold)                    % reset rounding
    end
  end
  