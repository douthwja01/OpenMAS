function r = mtimes(p,q)
%MTIMES       Polynomial multiplication  p * q
%
%p or q may be scalar (interval)
%
%Note that the result is always of type 'polynom', even if the result is constant.
%

% written  07/22/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
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
  
  if ~isa(p,'polynom')                            % p is scalar
      r = q;
      r.c = r.c * p;  
      setround(rndold)
      return
  end
  
  if ~isa(q,'polynom')                            % q is scalar
      r = p;
      r.c = r.c * q;  
      setround(rndold)
      return
  end
  
  if ( size(p.e,2)<=1 ) & ...                     % p univariate
          ( size(q.e,2)<=1 ) & ...                    % q univariate
          isequal(p.v,q.v)                            % variable of p and q identical
      r.e = p.e + q.e;
      if ~isa(p.c,'intval') & ~isa(q.c,'intval')    % non-interval data
          r.c = conv(p.c,q.c);
      else                                          % interval data
          rc = p.c.' * q.c;
          index = fliplr(toeplitz(p.e:-1:0,p.e:p.e+q.e));
          if isreal(rc)
              setround(-1)
              rcinf = flipud(full(sparse(index(:)+1,1,rc.inf(:)))).';
              setround(1)
              rcsup = flipud(full(sparse(index(:)+1,1,rc.sup(:)))).';
              r.c = intval(rcinf,rcsup,'infsup');      
          else
              rcmid = rc.mid(:);
              setround(-1)
              rcmidinf = flipud(full(sparse(index(:)+1,1,rcmid))).';
              setround(1)
              rcmidsup = flipud(full(sparse(index(:)+1,1,rcmid))).';
              rcrad = flipud(full(sparse(index(:)+1,1,rc.rad(:)))).';
              rcmid = rcmidinf + 0.5*(rcmidsup-rcmidinf);
              rcrad = abs(rcmid-rcmidinf) + rcrad;
              r.c = midrad(rcmid,rcrad); 
          end
      end
      r.v = p.v;
      r = class(r,'polynom');  
      setround(rndold)
      return
  end
  
  % result p*q multivariate
  np = length(p.c);                     % number of terms in p
  nq = length(q.c);                     % number of terms in q
  if ~isequal(p.v,q.v)                  % identical set of variables, both p and q must be multivariate
      [p.v,I,J] = joinvars(p.v,q.v);
      n = length(p.v);                    % number of variables
      pe = p.e;
      p.e = zeros(np,n);
      if size(pe,2)>1                     % p multivariate
          p.e(:,I) = pe;
      else                                % p univariate
          p.e(:,I) = (np-1:-1:0)';
          p.c = p.c.';
      end
      qe = q.e;
      q.e = zeros(nq,n);
      if size(qe,2)>1                     % q multivariate
          q.e(:,J) = qe;
      else                                % q univariate
          q.e(:,J) = (nq-1:-1:0)';
          q.c = q.c.';
      end
  else
      n = length(p.v);                    % number of variables
  end
  
  if np==1                              % only one term in p
      r.e = ones(nq,1)*p.e + q.e;
      r.c = p.c * q.c;
  elseif nq==1                          % only one term in q
      r.e = p.e + ones(np,1)*q.e;
      r.c = p.c * q.c;
  else
      pe = repmat(p.e,[1 1 nq]);
      qe = shiftdim(repmat(q.e',[1 1 np]),2); 
      r.e = reshape(shiftdim(pe + qe,2) , np*nq , n );
      r.c = reshape( q.c * p.c.' , np*nq , 1 );
  end
  r.v = p.v;
  
  r = normalize( class(r,'polynom') );
  
  setround(rndold)
