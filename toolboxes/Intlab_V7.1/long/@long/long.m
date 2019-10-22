function C = long(d)
%LONG         Long class constructor
%
%  C = long(d)
%
%  C is long scalar or column vector, array input d forced to be column vector
%
%Long representation
%
%  C = C.sign * sum( C.mantissa(i)*beta^(-i) ) * beta^C.exponent +/- C.error
%
%where C.error = C.error.mant * beta^C.error.exp and with
%  summation from 1 to precision (=size(C.mantissa,2)) with
%
%  1 <= precision <= INTLAB_LONG_PRECISION
%
%and integers
%
%  C.sign       in {-1,1}
%  C.mantissa   in 0 .. beta-1
%  C.exponent   representable integer ( -2^52+1 .. 2^52-1 )
%  C.error      nonnegative double, stored by C.error.mant and C.error.exp
%
%Computations can be executed with or w/o error term, see longinit. For
%  computational speed, comparison, min/max and absolute value refer to
%  midpoint.
%  To compare, for example, intervals (for computation with error term)
%  use inf(A)>sup(B) instead of A>B, and so forth.
%
%For control of working precision, see help longprecision. Base beta is a
%  power of 2 so that double precision floating point numbers are stored
%  without error.
%An example of long arithmetic with big cancellation is
%
%  x = long(-20);
%  y = long(1);  t = long(1);  i = 0;
%  while abs(t)>1e-20
%    i = i+1;
%    t = t*x/i;
%    y = y+t;
%  end
%  format long
%  Y = long2intval(y)
%
%producing
%
%  intval Y =
%    1.0e-008 *
%     0.20612_________
%
%The poor accuracy improves with more precision. After
%  longprecision(40);
%and the same statements as above we obtain
%
%  intval Y =
%    1.0e-008 *
%     0.20611536224378
%
%For more information try demolong .
%

% written  12/30/98     S.M. Rump
% modified 02/09/01     S.M. Rump  performance improvement
% modified 09/29/02     S.M. Rump  care for NaN components
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  
  if nargin==0
    C.sign = [];
    C.exponent = [];
    C.mantissa = [];
    C.error.mant = [];
    C.error.exp = [];
    C = class(C,'long');
    return
  end

  if isa(d,'long')
    C = d;
  else
    if isempty(d)
      C.sign = [];
      C.exponent = [];
      C.mantissa = [];
      C.error.mant = [];
      C.error.exp = [];
      C = class(C,'long');
      return
    end
    sized = size(d);
    n = prod(sized);
    if n~=sized(1)
      warning('input array for long forced to be column vector')
      d = d(:);
    end
    indexNaN = isnan(d);
    d(indexNaN) = [];
    [s e m] = splitdble(d);

    % get sign
    C.sign = s;

    % get exponent
    q = ceil(e/INTLAB_LONG_LOGBETA);
    C.exponent = q;

    % get mantissa digits
    C.mantissa = zeros(size(d,1),ceil(53/INTLAB_LONG_LOGBETA)+1);
    p = 0;
    while any(m)
      p = p+1;
      m = m*INTLAB_LONG_BETA;
      C.mantissa(:,p) = floor(m);
      m = m - C.mantissa(:,p);
    end

    % treat zero components
    index = ( d==0 );
    if any(index)
      C.sign(index) = 1;
      C.exponent(index) = -inf;
    end

    % adjust mantissa digits to exponent
    r = INTLAB_LONG_LOGBETA*q - e;
    index = ( r~=0 );
    if any(index)
      C.mantissa(index,:) = shiftright(C.mantissa(index,:),r(index));
    end

    % omit trailing zeros (improves performance)
    [m index] = max(C.mantissa(:,end:-1:1)~=0,[],2);
    index(all(C.mantissa'==0)) = size(C.mantissa,2);
    if min(index)~=1
      C.mantissa = C.mantissa(:,1:end-min(index)+1);
    end

    if any(indexNaN)
      Csign = C.sign;
      Cexponent = C.exponent;
      Cmantissa = C.mantissa;
      C.sign = zeros(n,1);
      C.exponent = zeros(n,1);
      C.mantissa = zeros(n,size(Cmantissa,2));
      C.sign(indexNaN) = NaN;
      C.exponent(indexNaN) = NaN;
      C.mantissa(indexNaN) = NaN;
      C.sign(~indexNaN) = Csign;
      C.exponent(~indexNaN) = Cexponent;
      C.mantissa(~indexNaN) = Cmantissa;
    end
    
    % set error
    C.error.mant = zeros(n,1);
    C.error.exp = zeros(n,1);

    C = class(C,'long');

  end
