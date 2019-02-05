function res = relerr(c,d)
%RELERR       Entrywise relative error
%             c  may be interval scalar, vector or matrix
%
%   res = relerr(c)
%
% if 0 not in c     rad/abs(mid)
% if 0 in c         rad
%
%
%For two input arguments,
%
%   res = relerr(c,d)
%
%with  res = max(relerr(c.inf,d.inf),relerr(c.sup,d.sup))
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  two arguments
% modified 08/07/02     S.M. Rump  abs to abss, definition changed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/08/05     S.M. Rump  NaN corrected
% modified 12/05/06     S.M. Rump  sparse output
% modified 01/03/07     S.M. Rump  redesign
% modified 05/22/07     S.M. Rump  Matlab bugs removed
% modified 10/06/07     S.M. Rump  inf/NaN correction (thanks to F. Bornemann)
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if nargin==1
    Index = isnan(c);
    if any(Index(:))
      %VVVV  c(Index) = 0;            % take care of sparse input
      s.type = '()'; s.subs = {Index}; c = subsasgn(c,s,0);
      %AAAA  Matlab bug fix
    end
  else
    Index = isnan(c) | isnan(d);
    if any(Index(:))
      %VVVV  c(Index) = 0;  d(Index) = 0;
      s.type = '()'; s.subs = {Index}; c = subsasgn(c,s,0); d = subsasgn(d,s,0);
      %AAAA  Matlab bug fix
    end
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==1
    if c.complex
      [m,n] = size(c.mid);
      if issparse(c.mid)
        res = sparse([],[],[],m,n);
      else
        res = zeros(m,n);
      end
      if isempty(c.mid)
        if rndold
          setround(rndold)
        end
        return
      end
    else
      [m,n] = size(c.inf);
      if issparse(c.inf)
        res = sparse([],[],[],m,n);
      else
        res = zeros(m,n);
      end
      if isempty(c.inf)
        if rndold
          setround(rndold)
        end
        return
      end
    end
    cmid = mid(c);
    crad = abs(rad(c));
    index = ( crad~=0 );
    if any(index(:))
      abscmid = abs(cmid(index));
      crad_ = crad(index);
      N = (crad_<abscmid).*abscmid + (crad_>=abscmid);
      res(index) = crad_ ./ N;
    end
    index = isinf(cmid) | isinf(crad);
    if any(index(:))
      res(index) = inf;
    end
    if any(Index(:))
      res(Index) = NaN;            % take care of sparse input
    end
  else
    res = max(relerr(inf_(c),inf_(d)),relerr(sup(c),sup(d)));
    if any(Index(:))
      res(Index) = NaN;            % take care of sparse input
    end
  end
  
  if rndold
    setround(rndold)
  end
