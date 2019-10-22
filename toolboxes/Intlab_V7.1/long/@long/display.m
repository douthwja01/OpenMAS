function display(C,digits)
%DISPLAY      Command window display for long
%
%This routine is called automatically when leaving the semicolon off
%an expression. The call
%
%    display(C)
%
%has the same effect. Calling display with a second parameter
%
%    display(C,digits)
%
%displays approximately digits decimals of the long number (or array) C.
%For digits=0, all figures up to the accuracy of C are printed
%  (careful: slow for exponents far away from zero because powers
%   of 10 have to computed exactly)
%Actual error is taken into account, only last 5 or 6 figures may be incorrect.
%
%!!! Careful, rather than output of intval, long output is *not* verified !!!
%!!! To produce a true display, use long2intval(C) !!!
%

%For internal use, display(C,-1) displays internal components of C,
%display(C,-2) produces Maple output
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 09/20/04     S.M. Rump  comment changed
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
  
  if isempty(inputname(1))
    name = 'ans';
  else
    name = inputname(1);
  end

  if isempty(C.sign)
    disp([ 'long ' name ' = '])
    disp('     []')
    return
  end

  if ( nargin>1 ) & ( digits==-1 )
    disp([ 'long ' name '.sign = ']), disp(C.sign)
    disp([ 'long ' name '.exponent = ']), disp(C.exponent)
    disp([ 'long ' name '.mantissa = ']), disp(C.mantissa)
    if INTLAB_LONG_ERROR
      disp([ 'long ' name '.error = ']), disp([ C.error.mant C.error.exp ])
    end
  end

  if ( nargin==2 ) & ( digits==-2 )    % Maple output, internal use only
    precC = length(C.mantissa);
    disp([ name ':=array([' ])
    s=[];
    for i=1:precC-1
      s = [ s sprintf('%10.0f',C.mantissa(i)) ',' ];
      if mod(i,6)==0, disp(s), s=[]; end
    end
    disp([ s sprintf('%10.0f',C.mantissa(end))  ' ]):' ])
    disp([ 's:=0: for i from ' int2str(precC) ...
           ' by -1 to 1 do s:=' name '[i]+s/2^25; od:' ])
    disp(['evalf(s,' int2str(longprecision) ');' ])  
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  loose = strcmp(get(0,'FormatSpacing'),'loose');
  disp([ 'long ' name ' = ' ])

  % take maximum first 80 bits
  p = min( size(C.mantissa,2) , ceil(80/INTLAB_LONG_LOGBETA)+1 );
  indexzero = all( C.mantissa(:,1:p)==0 , 2 );
  ws = warning;
  warning off
  logC = log10( C.mantissa(:,1:p) * (INTLAB_LONG_BETA.^(-1:-1:-p))' ) + ...
            C.exponent * log10(INTLAB_LONG_BETA);
  warning(ws)
  q = floor( logC );

  if ( nargin==1 ) | ( digits==-1 )

    for i=1:size(C.mantissa,1)
      if indexzero(i)
        mant = 0;
      else
        mant = C.sign(i)*10^(logC(i)-q(i)) ;
      end
      disp( [ sprintf('%16.12f',mant) ' * 10^' sprintf('%+05d',q(i)) ] )
    end

  else

    if digits==0         %%%%%%%%indexnonzero
      digits = longprecision-1;
    end
    n = size(C.mantissa,1);
    logbeta10 = floor( log10(INTLAB_LONG_BETA) );
    beta10 = 10^logbeta10;
    format = [ '%0' int2str(logbeta10) '.0f' ];
    Ci = long;

    for i=1:n
%VVVV Ci = C(i);
      Ci.sign = C.sign(i);
      Ci.exponent = C.exponent(i);
      Ci.mantissa = C.mantissa(i,:);
      Ci.error.mant = C.error.mant(i);
      Ci.error.exp = C.error.exp(i);
%AAAA Matlab V5.2 bug fix
      d = 1;
      drow = 1;
      s = [];
      first = 1;
      if q(i)>0
        prec = longprecision;
        longprecision( prec+15 );
        % k=15 is maximum k with beta * 10^k < 2^52
        Q = floor(q(i)/15);
        R = q(i) - 15*Q;
        for j=1:Q
          Ci = Ci / 1e15;
        end
        if R~=0
          Ci = Ci / 10^R;
        end
        longprecision(prec);
      elseif ( q(i)<0 ) & isfinite(q(i))
        prec = longprecision;
        longprecision( min( abs(q(i))+10 , prec ) );
        pow10 = long(10)^(-q(i));
        longprecision(prec);
        Ci = Ci * pow10;
      end
      if n>1
        disp([ name '(' int2str(i) ') =' ])
      end
      while 1
        if all( Ci.mantissa==0 )
          intpart = 0;
        else
          intpart = Ci.mantissa(1:Ci.exponent) * ...
                       INTLAB_LONG_BETA.^(Ci.exponent-1:-1:0)';
        end
        if first
          s = [ s sprintf('%3.0f',Ci.sign*intpart) '.' ];
          Ci.sign = 1;
          first = 0;
        else
          s = [ s sprintf(format,intpart) ];
        end
        Ci = Ci - intpart;
%       if ( size(Ci.mantissa,2)==1 ) | ( d>=digits )
        if ( Ci.error.mant*(INTLAB_LONG_BETA^Ci.error.exp)>=1 ) | ( d>=digits )
          break
        end
        d = d + logbeta10;
        Ci = Ci * beta10;
        if drow>=70
          disp([ s ' \' ])
          s = '    ';
          drow = logbeta10;
        else
          drow = drow+logbeta10;
        end
      end
      s = [ s ' * 10^' sprintf('%+05d',q(i)) ];
      disp(s)

    end

  end

  if loose, disp(' '); end
  
  if rndold
    setround(rndold)
  end
