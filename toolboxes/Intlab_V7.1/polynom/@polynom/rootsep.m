function d = rootsep(P,real)
%ROOTSEP      minimum distance between distinct roots of P (numerical)
%
%  d = rootsep(P,real)
%
%defined by   
%  d = min{ abs(a-b) :  P(a)=P(b)=0 and a~=b }
%Result is purely numerical and may afflicted with heavy rounding errors,
%  especially in the presence of multiple roots.
%If second parameter "real" is specified, then roots a,b are real.
%

% written  04/13/08     S.M. Rump
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if numvars(P)>1
    error('rootsep only for univariate polynomials')
  end
  
  e = 1e-30;
  if 1+e==1-e                   % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  if nargin<2                   % set defaults
    real = 0;
  end
  
  R = roots(P);
  if real
    R = R(imag(R)==0);
  end
  if ~isempty(R)
    n = degree(P);
    R = R(:)*ones(1,length(R));
    d = abs(R-R.');
    d(d==0) = inf;
    d = min(d(:));
  else
    d = inf;
  end
  
  if rndold
    setround(rndold)
  end
  