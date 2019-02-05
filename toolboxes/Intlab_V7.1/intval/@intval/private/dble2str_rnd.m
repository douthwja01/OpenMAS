function str = dble2str_rnd(c,commonexp,len,prec,expon,rnd)
%DBLE2STR_RND   internal routine for output routines infsup, midrad
%               with rigorous rounding
%

%Real number "c" is rounded to output format %len.precX with
%  X=e/f for expon=1/0, respectively.
%In case expon=0 (i.e. f-format), 0 <= c <= 1000 required (format short or long)
%  commonexp~=0 only in case expon=0
%For output format: 2 <= prec <= 15 with len large enough to carry all digits
%Output string "str" is according to format such that content of
%  string "str" is c*10^commenexp correctly rounded according to "rnd"
%

% written  12/30/98     S.M. Rump
% modified 01/19/04     S.M. Rump  rounding to nearest when calling sprintf
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/04/12     S.M. Rump  take care of sparse input
%

  % exception handling
  if isnan(c)
    str = [ blanks(len-3) 'NaN' ];
    return
  end
  if isinf(c)
    str = [ blanks(len-3) 'Inf' ];
    if c < 0
      str(end-4) = '-';
    end
    return
  end

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');
  INTLAB_INTVAL_ETA = realmin*eps;

  crnd = full(abs(c));
  signc = sign(c);
  round = signc*rnd;

  factor = INTLAB_INTVAL_POWER10.sup(1,commonexp+341);
  pprec = min(prec+4,15);
%VVVV  bug fix for Watcom underflow printf
  if crnd<realmin
    pprec = 15;
  end
%AAAA  end bug fix
  formatstr = [ '%' sprintf('%d',pprec+7) '.' sprintf('%d',pprec) 'e' ];

  while 1
    setround(0)
    str = sprintf(formatstr,crnd/factor);
    if isequal(str(end-2:end),'Inf')
      str = [ blanks(len-3) 'Inf' ];
      if c < 0
        str(end-4) = '-';
      end
      return
    end
    % make sure that exponent has three digits
    if ( lower(str(end-3))=='e' )
      str = [ str(2:end-2) '0' str(end-1:end) ];
    end

    Exp = (str(end-2:end)-'0')*[100;10;1];
    if isequal(str(end-3) , '-')
      Exp = -Exp;
    end
    if expon               % e-format
      digits = prec;
    else                   % f-format
      digits = prec + Exp;
    end

    % convert str with "digits" figures after decimal point back to double
    m = str([digits+2:-1:3 1]) - '0';
    offset = 9*(commonexp+Exp+341);
    mant = 9*(digits+1:-1:1);
    index = (m~=0).*(offset - mant + m) + (m==0);
    setround(-round);
    if round==1
      xrnd = sum( INTLAB_INTVAL_POWER10.inf(index) );
    else
      xrnd = sum( INTLAB_INTVAL_POWER10.sup(index) );
    end

    setround(round);
    if ( round*(abs(c)-xrnd) <= 0 )
      if expon     % e-format
        str = [ blanks(len-prec-7) str(1:digits+2) str(end-4:end) ];
        if signc == -1
          str(len-prec-7) = '-';
        end
      else                 % f-format
        str(2) = '';       % omit decimal point (digits+1 mantissa digits)
        mantissa = str(1:Exp+prec+1);
        if Exp < 0         % make sure length(mantissa) >= prec+1
          mantissa = [ char('0'+zeros(1,prec+1-length(mantissa))) mantissa ];
        end
        mantissa = [ mantissa(1:end-prec) '.' mantissa(end-prec+1:end) ];
        if signc == -1
          mantissa = [ '-' mantissa ];
        end
        str = [ blanks(len-length(mantissa)) mantissa ];
      end
      break
    end
    e = commonexp+Exp-digits;
    if expon & round==-1 & all( str(3:digits+2)=='0' )
      e = e-1;
    end
    delta = INTLAB_INTVAL_POWER10.sup(1,e+341);   % 10^e
    delta = max(delta,INTLAB_INTVAL_ETA);
%   setround(round);       % round is set to -round
    crnd = crnd + round*delta;
  end
  setround(0)
