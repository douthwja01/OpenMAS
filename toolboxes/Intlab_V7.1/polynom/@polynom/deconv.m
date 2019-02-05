function [q,r] = deconv(p1,p2)
%DECONV       Implements  p1 / p2  for univariate polynomials with remainder term
%
%   [q,r] = deconv(p1,p2)
%

% written  07/24/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(p1,'polynom')
    p1 = polynom(p1,'');
  end
  
  if ~isa(p2,'polynom')
    p2 = polynom(p2,'');
  end
  
  if ( size(p1.e,2)>1 ) | ( size(p2.e,2)>1 )
    error('divison only for univariate polynomials')
  end
  
  % both p1 and p2 univariate polynomials
  if ~isequal(p1.v,p2.v) & ~isempty(p1.v) & ~isempty(p2.v)
    error('division only for univariate polynomials depending on the same variable')
  end
  
  isint = isa(p1.c,'intval') | isa(p2.c,'intval');
  
  if p2.e==0                           % p2 is constant
    q = p1;
    q.c = p1.c/p2.c;
    if nargout==2
      if isint
        r = polynom(intval(0),p1.v);
      else
        r = polynom(0,p1.v);
      end
    end  
    if rndold
      setround(rndold)
    end
    return
  end
  
  if p1.e<p2.e                          % deg(p1) < deg(p2)
    q.e = 0;
    if isint
      q.c = intval(0);
    else
      q.c = 0;
    end
    q.v = p1.v;
    q = class(q,'polynom');
    if nargout==2
      r = p1;
    end  
    if rndold
      setround(rndold)
    end
    return
  end    
  
  q.e = max(p1.e-p2.e,0);
  if isint
    q.c = intval(zeros(1,q.e+1));
    p1.c = intval(p1.c);
    for i=1:(q.e+1)
      q.c(i) = p1.c(i)/p2.c(1);
      p1.c(i:(i+p2.e)) = p1.c(i:(i+p2.e)) - q.c(i)*p2.c;
    end
    if nargout==2
      np1 = length(p1.c);
      r.e = p2.e-1;
      r.c = p1.c((np1-p2.e+1):np1);
      r.v = p1.v;
      r = class(r,'polynom');
    end
  else
    if nargout==2
      [q.c,rc] = deconv(p1.c,p2.c);
      r.e = min(p2.e-1,length(rc)-1);
      r.c = rc(end-r.e:end);  
      r.v = p1.v;
      r = class(r,'polynom');
    else
      q.c = deconv(p1.c,p2.c);
    end
  end
  q.v = p1.v;
  
  q = class(q,'polynom');
  
  if rndold
    setround(rndold)
  end
