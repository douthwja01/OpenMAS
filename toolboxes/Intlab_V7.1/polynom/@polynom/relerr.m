function res = relerr(p,q)
%RELERR       Maximum relative error between coefficients of polynomials p and q
%
%   res = relerr(p,q)
%

% written  09/02/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/05/06     S.M. Rump  one input argument
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  if nargin==1
    res = relerr(inf_(p),sup(p));
    return
  end

  p = polynom(p);
  q = polynom(q);

  np = length(p.c);
  nq = length(q.c);
  [z,I,J] = joinvars(p.v,q.v);

  % recompute expononents with respect to joint set of variables
  pe = zeros(np,length(z));
  if size(p.e,2)>1                  % p multivariate
    pe(:,I) = p.e;
  else                              % p univariate
    pe(:,I) = (np-1:-1:0)';
  end
  qe = zeros(nq,length(z));
  if size(q.e,2)>1                  % q multivariate
    qe(:,J) = q.e;
  else                              % q univariate
    qe(:,J) = (nq-1:-1:0)';
  end

  % p.c and q.c coefficients w.r.t. pe and qe, respectively.
  % pe and qe have same number of columns (corresponding to variables),
  %   but possibly different number of rows:
  % make sure pe and qe are identical

  pmax = max(pe,[],1);
  qmax = max(qe,[],1);
  emax = max(pmax,qmax)+1;
  ip = exp2index(pe,emax)+1;     % scalar indices of pe w.r.t. emax
  iq = exp2index(qe,emax)+1;     % scalar indices of qe w.r.t. emax

  % put coefficients p.c and q.c into linearized array
  n = max(max(ip),max(iq));
  pc = typeadj( zeros(n,1) , typeof(p.c) );
  pc(ip) = p.c;
  qc = typeadj( zeros(n,1) , typeof(q.c) );
  qc(iq) = q.c;

  res = max( relerr(pc,qc) );
  
  if rndold
    setround(rndold)
  end
