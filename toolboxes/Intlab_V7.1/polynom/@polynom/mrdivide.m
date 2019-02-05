function r = mrdivide(p,q)
%MRDIVIDE     Implements  p / q  for univariate polynomials
%

% written  07/24/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if ~isa(q,'polynom')                % p is constant
      r = p;
      r.c = p.c/q;
      return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(p,'polynom')
      if iscell(q.v)
          error('division only for univariate polynomials')
      end
      if q.e==0                         % q is univariate constant
          r = q;
          r.c = p/q.c;  
          if rndold
            setround(rndold)
          end
          return
      end
      r.e = 0;
      r.c = typeadj(0,typeof(q.c));
      r.v = q.v;
      r = class(r,'polynom');  
      if rndold
        setround(rndold)
      end
      return
    end
  
  if ( size(p.e,2)>1 ) | ( size(q.e,2)>1 )
      error('divison only for univariate polynomials')
  end
  
  % both p and q univariate polynomials
  if ~isequal(p.v,q.v) & ~isempty(p.v) & ~isempty(q.v)
      error('division only for univariate polynomials depending on the same variable')
  end
  
  if q.e==0                           % q is constant
      r = p;
      r.c = p.c/q.c;  
      if rndold
        setround(rndold)
      end
      return
  end
  
  isint = isa(p.c,'intval') | isa(q.c,'intval');
  
  if p.e<q.e                          % degree(p) < degree(q)
      r.e = 0;
      if isint
          r.c = intval(0);
      else
          r.c = 0;
      end
      r.v = q.v;
      r = class(r,'polynom');  
      if rndold
        setround(rndold)
      end
      return
  end    
  
  r.e = max(p.e-q.e,0);
  if isint
      r.c = intval(zeros(1,r.e+1));
      p.c = intval(p.c);
      for i=1:(r.e+1)
          r.c(i) = p.c(i)/q.c(1);
          p.c(i:(i+q.e)) = p.c(i:(i+q.e)) - r.c(i)*q.c;
      end
  else
      r.c = deconv(p.c,q.c);
  end
  r.v = p.v;
  
  r = class(r,'polynom');
  
  if rndold
    setround(rndold)
  end
