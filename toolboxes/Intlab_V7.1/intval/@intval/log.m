function y = log(x)
%LOG          Implements  log(x)  for intervals
%
%   y = log(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, following N.C. Boersken:
%                                  Komplexe Kreis-Standardfunktionen,
%                                  Freiburger Intervallberichte 78/2,
%                                  NaN input, sparse input, log(0),
%                                  major revision, improved accuracy
% modified 12/06/99                branch cut with warning
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and some
%                                     improvements, tocmplx replaced by cintval,
%                                     extreme values for approximate part
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 10/23/07     S.M. Rump  complex numbers
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/20/08     S.M. Rump  check for zero omitted, improved performance
% modified 10/18/08     S.M. Rump  StdFctsException ignore/NaN
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  if x.complex
    if issparse(x.mid)
      x.mid = full(x.mid);
      x.rad = full(x.rad);
    end
    y = x;

%   y.mid = log(x.mid);
%   y.rad = - log( 1 - x.rad./abs(x.mid) );

    INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

    xmidre = real(x.mid);
    xmidim = imag(x.mid);
    Xmidim = intval(xmidim);
    Mim = atan(Xmidim./xmidre);      % -pi/2 <= Mim <= pi/2

    % special treatment of imaginary axis
    index = ( xmidre==0 );
    if any(index(:))
      Mim.inf(index) = INTLAB_STDFCTS_PI.PI2INF;          % pi/2
      Mim.sup(index) = INTLAB_STDFCTS_PI.PI2SUP;
      indexneg = index & ( xmidim<0 );
      if any(indexneg(:))
        Mim.inf(indexneg) = -INTLAB_STDFCTS_PI.PI2SUP;    % -pi/2
        Mim.sup(indexneg) = -INTLAB_STDFCTS_PI.PI2INF;
      end
    end

    index = ( xmidre < 0 );
    if any(index(:))
      % correct to atan2:  -pi <= Mim <= pi
      Pi = infsup( INTLAB_STDFCTS_PI.PIINF, ...
                   INTLAB_STDFCTS_PI.PISUP );
      signxmidim = sign(xmidim(index));
      corr = Pi.*(signxmidim+(signxmidim==0));
      Mim.inf(index) = Mim.inf(index) + corr.inf;
      Mim.sup(index) = Mim.sup(index) + corr.sup;
    end

    % x.mid = r*exp(j*phi)  ==>  log(x.mid) = log(r) + j*phi
    R = abs(x.mid);
    Mreinf = log_rnd(R,-1);
    Mresup = log_rnd(R,1);

    setround(1)
    mre = Mreinf + 0.5*(Mresup-Mreinf);
    mim = Mim.inf + 0.5*(Mim.sup-Mim.inf);
    y.mid = mre + j*mim;
    mrad = abs( mre-Mreinf + j*(mim-Mim.inf) );
    % log(R)  in  mre + j*mim +/- mrad

    wng = warning;
    warning off
    % x.rad < R ,  otherwise zero interval
    y.rad = log_rnd( -( x.rad./R - 1 ) , -1);
    warning(wng);
    y.rad = -y.rad + mrad;
    setround(0)

    % special treatment of branch cut
    index0 = ( R <= x.rad );               % zero interval
    index = ( ~index0 ) & ...
            (   ( ( xmidre<0 ) & ( xmidim>=0 ) & ( x.rad>xmidim ) ) ...
              | ( ( xmidre<0 ) & ( xmidim<0 ) & ( x.rad>=-xmidim ) ) ...
            );
    if any(index(:))
      warning('Complex Log: Input interval intersects with branch cut')
    end

    if any(index0(:))
      y.mid(index0) = complex(NaN,NaN);
      y.rad(index0) = NaN;
    end
  
    if rndold
      setround(rndold)
    end
    
    return
  end

  if issparse(x.inf)
    x.inf = full(x.inf);
    x.sup = full(x.sup);
  end
  % input x real and full
  % real range of definition:  [0,inf]
  INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
  index = ( x.inf<0 );                  % (partially) exceptional indices
  if any(index(:))                      % handle input out-of-range
    if INTLAB_STDFCTS_EXCPTN<=1  % out-of-range input handled as complex
      if INTLAB_STDFCTS_EXCPTN==1
        warning('LOG: Real interval input out of range changed to be complex')
      end
      y = x;
      %VVVV  y(index) = log(cintval(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,log(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~index;
      if any(index(:))
        %VVVV  y(index) = log(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,log(subsref(x,s)));
        %AAAA  Matlab bug fix
      end
      if rndold
        setround(rndold)
      end
      return
    end
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',1);
    if INTLAB_STDFCTS_EXCPTN==3    % ignore input out of range (ignore-mode)
      x.inf(index) = realmin*eps;;        % completely exceptional indices treated below
      indexneg = index & ( x.sup<0);      % completely exceptional indices
    end
  else
    indexneg = [];                        % make sure indexneg is not undefined
  end

  % input x real and full
  y = x;
  wng = warning;
  warning off
  
  % treat non-exceptional cases
  y.inf = log_rnd(x.inf,-1);
  y.sup = log_rnd(x.sup,1);

  if INTLAB_STDFCTS_EXCPTN==3      % ignore input out of range (ignore-mode)
    if ~isempty(find(indexneg))           % completely exceptional arguments to NaN
      y.inf(indexneg) = NaN;
      y.sup(indexneg) = NaN;
    end
  else                                    % any input out of range to NaN (NaN-mode)
    if ~isempty(find(index))              % exceptional arguments to NaN
      y.inf(index) = NaN;
      y.sup(index) = NaN;
    end
  end

  warning(wng)
  setround(rndold)
