function y = asinh(x)
%ASINH        Implements  asinh(x)  for intervals
%
%   y = asinh(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, large arguments improved,
%                                  sparse input, small arguments improved,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 01/20/03     S.M. Rump  Matlab sqrt fixed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  extreme values for approximate part
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,asinh(full(sx)),m,n);
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
    y = log( x + sqrt(sqr(x)+1) ); 
    if rndold
      setround(rndold)
    end
    return
  end
  
  y = x;

  xinf = x.inf(:);
  xsup = x.sup(:);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup(:)<=0 );
  len2 = sum(IndexSupNeg);

  Y = asinhpos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ~IndexInfPos;
  len1 = sum(IndexInfNeg);
  IndexSupPos = ~IndexSupNeg;
  len2 = sum(IndexSupPos);

  Y = asinhpos( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : len1+len2 );

  setround(rndold)


function y = asinhpos(x,rnd)
% rigorous asinh(x) for nonnegative double vector x with
% rounding corresponding to rnd

  y = x;

  % huge arguments: use proper scaling
  index = ( x>1e10 );
  if any(index)            % x > 1e10
    setround(rnd)
    E = 2*x(index) + eps*(rnd+1);
    Y = log_rnd( E , rnd );
    y(index) = Y;
  end
  Index = ~index;

  % large arguments: use asinh(x) = log( 1 + ( 2*x*sqrt(x^2+1)+2*x^2 ))/2
  index = Index & ( x>0.35 );
  if any(index)            % 0.35 < x <= 1e10
    X = x(index);
    setround(rnd)
    XX = X.*X;
    E = 1 + ( 2*X.*sqrt_rnd( XX+1 , rnd ) + 2*XX );
    y(index) = log_rnd(E,rnd) / 2;

  end
  Index = Index & ( ~index );

  % small arguments: asinh(x) = log( 1 + 2*sqrt(x^4+x^2) + 2*x^2 )/2
  index = Index & ( x>1e-9);
  if any(index)            % 1e-9 <= x <= 0.35 ,  E must be <=1
    X = x(index);
    XX = X.*X;
    setround(rnd)
    E = 2*sqrt_rnd( XX.*XX+XX , rnd ) + 2*XX;
    y(index) = log_1(E,rnd) / 2;       % 0 <= E < 1
  end

  % very small arguments
  index = Index & ( ~index );          % 0 <= x <= 1e-9
  if any(index)
    X = x(index);
    setround(rnd)
    if rnd==-1                         % 0 <= err <= x^3/6 < 1.7e-19*x
      y(index) = X - 1.7e-19*X;
    else
      y(index) = X;
    end
  end
