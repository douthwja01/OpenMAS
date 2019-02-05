function [ X , Xintinf , Xintsup ] = structlss(As,b)
%STRUCTLSS    Dense linear system solver for structured matrices
%
%   [ X , Xintinf , Xintsup ] = structlss(As,b)
%
%input   As       Structured matrix, see below
%        b        right hand side, must be one column
%output  X        Inclusion of solution of structured linear system
%        Xintinf  (optional) inner inclusion of solution of structured
%        Xintsup    linear system
%
%Input As characterizes a structured matrix, see function "structure".
%Input must be real; for complex systems use augmented linear system.
%For structured right hand side use use augmented linear system.
%
%For example, a linear system with Hilbert matrix afflicted with absolute
%  tolerances 1e-10 and r.h.s. ones(n,1) can be solved as follows:
%
%    n = 5;
%    A = midrad( hilb(n) , 1e-10 );
%    As = structure(A,'symmetric');
%    X = structlss( As , ones(n,1) )
%
%For other structures, see function "structure".
%
%Linear systems may be treated subject to different structures producing
%  quite different results. Try, for example,
%
%    n = 4;  e = 1e-3; intvalinit('displayinfsup');
%    A = midrad( toeplitz([0 1+e e 1+e]),1e-4); b = A.mid*ones(n,1);
%    X1 = verifylss(A,b);
%    X2 = structlss(structure(A,'symmetric'),b);
%    X3 = structlss(structure(A,'symmetricToeplitz'),b);
%    X4 = structlss(structure(A,'generalToeplitz'),b);
%    X5 = structlss(structure(A,'persymmetric'),b);
%    X6 = structlss(structure(A,'circulant'),b);
%    res = [ X1 X2 X3 X4 X5 X6 ], rad(res)
%
%Careful with output of inner inclusions. To see correct numbers use
%
%    intval(Xintinf), intval(Xintsup)
%
%and look at the inner bounds. If Xintinf_i>Xintsup_i, then no inner inclusion
%  could be computed for i-th component, possibly due to wide input intervals
%  and/or ill-conditioned matrix.
%
%
%The algorithm generalizes
%  C. Jansson: Interval Linear Systems with Symmetric Matrices,
%    Skew-Symmetric Matrices and Dependencies in the Right Hand Side,
%    Computing 46, 265-274 (1991)
%to general linear matrix structures, see also
%  S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic.
%    Acta Numerica, 19:287-449, 2010. 
%

% written  07/26/99     S.M. Rump
%

  if ~isstruct(As)
    error('input matrix must be structured')
  end

  [n,k] = size(As.Phi);
  n = round(sqrt(n));

  A = reshape(As.Phi*As.p,n,n);
  midA = mid(A);
  midb = mid(b);

  % preconditioner: approximate inverse
  R = inv( midA ) ;

  % apply one residual iteration to approximate solution to
  %   ensure backward stability (Skeel)
  xs = R * midb ;
  xs = xs + R*(midb - midA*xs);

  maxmem = 4.1e6;      % split computation into matrices of limited size
                       % here maximum 4.1 Mbyte, enough for 100x100 symmetric

  % interval iteration
  % Z = R * (b - A*xs) ;
  %   = R*b - ( kron(xs',R)*Phi ) * As.p
  %   = R*b - ( R*( kron(xs',eye(n) ) * Phi ) * As.p

  [i,j,s] = find(As.Phi);
  q = floor((i-1)/n);
  if 8*n*k<=maxmem
    R_Phi_x = intval(R) * sparse(i-n*q,j,s.*xs(q+1),n,k) ;
    Z = R*intval(b) - R_Phi_x * As.p;
  else
    Z = intval(zeros(n,1));
    step = ceil(maxmem/(8*k));    % number of columns fitting into maxmem
    for ii=1:step:n
      iiend = min(ii+step-1,n);
      Z(ii:iiend) = ( intval(R(ii:iiend,:)) * sparse(i-n*q,j,s.*xs(q+1),n,k) ) ...
                       * As.p;
    end
    Z = R*intval(b) - Z;
  end

  C = speye(n) - R*A;
  Y = Z;
  E = 0.1*rad(Y)*hull(-1,1) + midrad(0,10*realmin);
  i = 0; imax = 15; ready = 0;
  while ~ready & i<imax
    i = i+1;
    X = Y + E;
    delta = C * X;
    Y = Z + delta;
    ready = all(all(in0(Y,X)));
  end

  if ready
    X = xs + Y;
    if nargout==3       % Calculate inner inclusion
      [Rbinf,Rbsup] = intprod(R,mid(b),rad(b));

      [pmid,prad] = intmidrad(As.p.inf,As.p.sup);

      setround(-1)
      CCsup = kron(xs',R)*As.Phi;   % inner inclusion of kron(xs',R)*Phi
      setround(1)
      CCinf = kron(xs',R)*As.Phi;
      [CCmid,CCrad] = intmidrad(CCinf,CCsup);

      setround(1)
      CCpmidinf = CCmid*pmid;
      setround(-1)
      CCpmidsup = CCmid*pmid;

      CCprad = CCrad*abs(pmid) + abs(CCmid)*prad + CCrad*(-prad);
      CCpsup = CCpmidsup + CCprad;
      setround(1)
      CCpinf = CCpmidinf - CCprad;

      Zsup = Rbsup - CCpinf;
      setround(1)
      Zinf = Rbinf - CCpsup;

      setround(1)
      Xintinf = xs + Zinf + delta.sup;
      setround(-1)
      Xintsup = xs + Zsup + delta.inf;

    end
  else
    X = repmat(NaN,size(X));
    Xintinf = +inf*ones(n,1);
    Xintsup = -inf*ones(n,1);
  end


function [xmid,xrad] = intmidrad(xinf,xsup)
% produces xinf <= xmid-xrad <= xmid+xrad <= xsup ,
%   i.e.  in( <xmid,xrad> , [xinf,xsup] ) = 1

  setround(-1)
  xmid = xinf + 0.5*(xsup-xinf);
  xrad = xmid - xinf;


function [yinf,ysup] = intprod(R,xmid,xrad)
% definitely for all i there are x1,x2 in xmid+/-xrad with
%   (R*x1)_i <= yinf_i  and  ysup_i <= (R*x2)_i
% i.e.
%   in( [yinf,ysup]  ,  R [*] <xmid,xrad> ) = 1

  setround(-1)
  yrad = abs(R)*xrad;
  ysup = R*xmid + yrad;
  setround(1)
  yinf = R*xmid - yrad;
