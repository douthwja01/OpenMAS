function c = infsup(a,b)
%INFSUP       Initialization of interval by infimum and supremum
%  computed such that [a,b] is enclosed in interval  c
%
%   c = infsup(a,b)
%

% written  10/16/98     S.M. Rump
% modified 12/06/98     S.M. Rump  NaN test added
% modified 06/24/99     S.M. Rump  check sparsity, multi-dimensional arrays
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/07/04     S.M. Rump  very large linear index
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/02/05     S.M. Rump  large indices
% modified 09/10/07     S.M. Rump  performance, huge arrays
% modified 10/14/08     S.M. Rump  handle [inf,inf]
% modified 02/28/12     S.M. Rump  check input is numeric (thanks to Joshua Ross, Adelaide)
% modified 08/26/12     S.M. Rump  global variables removed
%

  if ~isnumeric(a) | ~isnumeric(b)
    error('both bounds for infsup must be numeric')
  end

  if issparse(a)~=issparse(b)
    error('midpoint and radius must both be or none of them be sparse')
  end

  if ~isreal(a) | ~isreal(b)        % complex interval
    index = any( ( real(a)>real(b) ) | ( imag(a)>imag(b) ) );
    if any(index(:))
      error('improper intervals in call of infsup')
    end
    INTLAB_INTVAL_CINFSUPASGN = getappdata(0,'INTLAB_INTVAL_CINFSUPASGN');
    if INTLAB_INTVAL_CINFSUPASGN
      warning('complex interval defined by infsup(zinf,zsup) will be stored in mid/rad causing overestimation')
    end
    rndold = getround;
    setround(1)
    cmid = a + 0.5*(b - a);
    crad = max( abs(cmid-a) , abs(cmid-b) ) ;
    setround(rndold)
    index = isinf(a) | isinf(b);    % does not work for huge arrays
    if any(index(:))                % take care of [inf,inf]
      cmid(index) = 0;
      crad(index) = inf;
    end
    c = intval(cmid,crad,'midrad');
  else                              % real interval (inf/sup representation)
    index = any( a>b );
    if any(index(:))                % take care of very large linear index
      error('improper intervals in call of infsup')
    end
    c = intval(a,b,'infsup');
  end

  index = isnan(a) | isnan(b);
  nanindex = any(index);
  if any(nanindex(:))               % take care of very large linear index
    c(index) = NaN;
  end
