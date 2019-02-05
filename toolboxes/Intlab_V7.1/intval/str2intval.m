function [x,err] = str2intval(str)
%STR2INTVAL   Rigorous conversion string to intval
%
%Rigorous conversion of string into intval
%
%  x = str2intval(str)
%
%Input constant syntax as Matlab constants, e.g. +3.75e-2+4.65e+2i
%  exponent character 'e' or 'E'
%  imaginary unit character 'i' or 'j'
%
%Input constants may be afflicted with tolerances, like 5.4-3.14__j
%  where the interval is defined by adding 1 to and subtracting 1 from
%  the last displayed figure before the first "_". 
%Note that complex intervals are always stored in mid-rad representation, so 
%  X = intval('5.4-3.14__j') produces roughly a circle with midpoint 5.4-3.14
%  and radius 0.01
%Inputs specified without tolerance are treated as if infinitely many
%  zeros follow. Input must be finite.
%
%To avoid error messages use
%
%  [x,err] = str2intval(str)
%
%where  err =  0   normal end, no error detected
%              1   improper interval (lower bound greater than upper)
%              2   other errors
%
%In case of error, result x is NaN for safety.
%
%Input string may be one- or two-dimensional
%  output x is always column vector
%  careful with complex data when you need a row vector: use x.'
%  to obtain a different size, use reshape
%

% written  10/10/98     S.M. Rump
% modified 09/01/00     S.M. Rump  intval('1e-1000') corrected, rounding unchanged after use
% modified 09/29/02     S.M. Rump  fix due to different behaviour of logical of Matlab 6.5
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 01/01/05     S.M. Rump  error parameter added (thanks to Arnold)
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/17/08     S.M. Rump  take care of other languages
% modified 10/27/08     S.M. Rump  huge numbers, inf, NaN
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  try                                       % capture possible errors
    
    lasterr('');                            % reset lasterr
    err = 0;

    [M N] = size(str);
    t = str';                               % make sure string is
    t(size(t,1)+1,:) = blanks(size(t,2));   % one row of numbers
    str = t(:)';

    firstinfsup = find( str=='[' );         % find input intervals [ , ]
    if ~isempty(firstinfsup)
      lastinfsup = find( str==']' );
      index = zeros(1,length(str));
      index(firstinfsup) = 1;
      index(lastinfsup) = -1;
      index = cumsum(index);
      index(lastinfsup) = 1;                % index==1 for intervals [ , ]
      indexspace = ( isspace(str) & index ) ;
      str( indexspace ) = '';              % delete blanks within intervals
      index( indexspace ) = '';
      sepinfsup = find( ( str==',' ) & index );  % index of separating comma
      lastinfsup = find( str==']' );       % index of ']'

      if ( length(firstinfsup)~=length(lastinfsup) ) | ...
          ( length(firstinfsup)~=length(sepinfsup) ) | ...
          any( sepinfsup >= lastinfsup )
        error('invalid input for str2intval')
      end

      str(sepinfsup) = ' ';                % replace ',' separator and '>' for
      str(lastinfsup) = ' ';               % intervals <...,...> by blank
    end

    firstmidrad = find( str=='<' );        % find input intervals < , >
    if ~isempty(firstmidrad)
      lastmidrad = find( str=='>' );
      index = zeros(1,length(str));
      index(firstmidrad) = 1;
      index(lastmidrad) = -1;
      index = cumsum(index);
      index(lastmidrad) = 1;               % index==1 for intervals < , >
      indexspace = ( isspace(str) & index ) ;
      str( indexspace ) = '';              % delete blanks within intervals
      index( indexspace ) = '';
      sepmidrad = find( ( str==',' ) & index );  % index of separating comma
      lastmidrad = find( str=='>' );       % index of '>'

      if ( length(firstmidrad)~=length(lastmidrad) ) | ...
          ( length(firstmidrad)~=length(sepmidrad) ) | ...
          any( sepmidrad >= lastmidrad )
        error('invalid input for str2intval')
      end

      str(sepmidrad) = ' ';                % replace ',' separator and ']' for
      str(lastmidrad) = ' ';               % intervals [...,...] by blank
    end

    zero = [0 find(isspace(str))];         % index vector of blanks in string
    index = find(diff(zero)~=1);           % difference of index vector of blanks

    I=cumsum(diff(zero));
    first = zero(index)+1;
    last = zero(index+1)-1;

    if ~isempty(firstinfsup)
      indexinfsup = find( str(first)=='[' );
      str(first(indexinfsup)) = ' ';
      first(indexinfsup) = first(indexinfsup) + 1;
    end

    if ~isempty(firstmidrad)
      indexmidrad = find( str(first)=='<' );
      str(first(indexmidrad)) = ' ';
      first(indexmidrad) = first(indexmidrad) + 1;
    end

    x = str2IV(str,first,last);

    if ~isempty(firstinfsup)
      x(indexinfsup) = infsup(x(indexinfsup).inf,x(indexinfsup+1).sup);
    end

    if ~isempty(firstmidrad)
      setround(1)
      R = abs( x(indexmidrad).sup-x(indexmidrad).inf ) + ...
        abs( x(indexmidrad+1).sup );
      setround(0)
      x(indexmidrad) = midrad(x(indexmidrad).inf,R);
    end

    index = [];
    if ~isempty(firstinfsup)
      index = [ index indexinfsup+1 ];
    end
    if ~isempty(firstmidrad)
      index = [ index indexmidrad+1 ];
    end

    if ~isempty(index)
      x(index) = [];
    end

  catch                                  % check for errors

    if ~isempty(lasterr)
      if nargout<2
        error(lasterr)
      else
        s = lasterr;
        if findstr(s,'improper')
          err = 1;
        else
          err = 2;
        end
        x = NaN;
      end
    end

  end
  
  setround(rndold)

  

function x = str2IV(str,first,last)
%local function; "str" is one row and contains only numbers, no intervals
%  starting at str(first) and ending at str(last)

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');

  if isempty(first)
    error('Blank or empty input for str2intval')
  end

  n = length(last);                      % number of inputs
  sign = ones(n,1);                      % calculate vector of signs
  s = str(first);
  v = s=='+';
  str(first(v)) = ' ';
  first(v) = first(v) + 1;
  v = s=='-';
  str(first(v)) = ' ';
  sign(v) = -1;
  first(v) = first(v) + 1;
  % sign(i)*str2num(str(first(i):last(i))) = inputnumber(i)  for i=1:n

  s = str(last);
  % index of complex inputs
  cindex = ( s == 'i' ) | ( s == 'j' ) ;

  if any(cindex)       % complex input

    xmid(n,1) = j;     % sqrt(-1)
    xrad(n,1) = 0;
    rindex = ~cindex;
    % sign(i)*str2num(str(first(i):last(i))) = real inputnumber(i)
    %   for i in rindex

    if any(rindex)

      % treat real input
      realstr = extractstr(str,first(rindex),last(rindex));

      [ mantissa , exponent , tol ] = ...
         strmantexp2dble(realstr,first(rindex),last(rindex));

      xinf = dec2dble(sign(rindex),mantissa,exponent,-1);
      xsup = dec2dble(sign(rindex),mantissa,exponent,+1);
      setround(1)
      xmid(rindex) = xinf + 0.5*(xsup-xinf);
      xrad(rindex) = xmid(rindex) - xinf;
      xrad(rindex) = xrad(rindex) + INTLAB_INTVAL_POWER10.sup(1,tol)';
      setround(0)

    end

    % treat complex input
    firstcindex = first(cindex);
    lastcindex = last(cindex) - 1;
    signcindex = sign(cindex);
    % sign(i)*str2num(str(first(i):last(i))) = complex inputnumber(i)
    %   without imaginary unit character   for i in cindex

    cmplxstr = extractstr(str,firstcindex,lastcindex);
    xmidcmplx(sum(cindex),1) = j;  % sqrt(-1)
    xradcmplx(sum(cindex),1) = 0;

    % split complex string (cmplxstr,first(cindex),last(cindex)) into
    %   real and imaginary part
    input = str2matrix(cmplxstr,firstcindex,lastcindex,0)'; % complex inputs right bound
    index = ( input == '+' ) | ( input == '-' ) ;
    sumindex = sum(index,1);
    if any(sumindex>3)
      error('Invalid input for str2intval')
    end

    %%% leading sign already stored in "sign"
    iImagonly = ( sumindex==0 );     % no sign: input like 3e1i
                                     % definitely only imaginary part

    index1 = ( sumindex==1 );        % one sign: input either 1e+4i or 3+4i
    if any(index1)
      str1 = extractstr(cmplxstr,firstcindex(index1),lastcindex(index1));
      index = find( ( str1 == '+' ) | ( str1 == '-' ) );
      str1_1 = str1(index-1);
      findex1 = find(index1);
      iImagonly( findex1( find( (str1_1=='e') | (str1_1=='E') ) ) ) = 1;
    end
    % iImagonly: numbers with only imaginary part
    % no or one sign:
    % 4i or 3e+4i

    if any(iImagonly)

      strI = extractstr(cmplxstr,firstcindex(iImagonly),lastcindex(iImagonly));
      [ mantissa , exponent , tol ] = ...
         strmantexp2dble(strI,firstcindex(iImagonly),lastcindex(iImagonly));

      xinf = dec2dble(signcindex(iImagonly),mantissa,exponent,-1);
      xsup = dec2dble(signcindex(iImagonly),mantissa,exponent,+1);
      setround(1)
      xmidcmplx(iImagonly) = xinf + 0.5*(xsup-xinf);
      xradcmplx(iImagonly) = xmidcmplx(iImagonly) - xinf ...
                            + INTLAB_INTVAL_POWER10.sup(1,tol)';
      setround(0)
      xmidcmplx(iImagonly) = xmidcmplx(iImagonly) * j;   % make complex

    end

    iRealimag = ~ iImagonly;   % iRealimag: numbers with real and imaginary part
                               % one to three signs:
                               % 3e1+4i or 3e1+4e+2i or 3e+1+4e+2i
    if any(iRealimag)

      strRI = extractstr(cmplxstr,firstcindex(iRealimag),lastcindex(iRealimag));
      index = find( ( strRI == '+' ) | ( strRI == '-' ) );
      strRI_1 = strRI(index-1);
      % index of sign between real and imaginary part
      signindex = index( ( (strRI_1>='0') & (strRI_1<='9') ) ...
                           | (strRI_1=='.') | (strRI_1=='_') );

      % real part
      strre = extractstr(str,firstcindex(iRealimag),signindex-1);
      [ mantre , expre , tolre ] = ...
         strmantexp2dble(strre,firstcindex(iRealimag),signindex-1);

      % imaginary part
      strim = extractstr(str,signindex+1,lastcindex(iRealimag));
      [ mantim , expim , tolim ] = ...
         strmantexp2dble(strim,signindex+1,lastcindex(iRealimag));

      signim = 2*( strRI(signindex) == '+' ) - 1;
      xinf = dec2dble(signcindex(iRealimag),mantre,expre,-1) + ...
             dec2dble(signim,mantim,expim,-1)*j;
      xsup = dec2dble(signcindex(iRealimag),mantre,expre,+1) + ...
             dec2dble(signim,mantim,expim,+1)*j;
      setround(1)
      xmidcmplx(iRealimag) = xinf + 0.5*(xsup-xinf);
      tol = max(tolre,tolim);
      xradcmplx(iRealimag) = abs( xmidcmplx(iRealimag) - xinf ) ...
                              + INTLAB_INTVAL_POWER10.sup(1,tol)';
      setround(0)

    end

    % collect complex inputs
    xmid(cindex) = xmidcmplx;
    xrad(cindex) = xradcmplx;
    x = cintval(xmid,xrad);     % make sure output is complex

  else               % only real input

    [ mantissa , exponent , tol ] = strmantexp2dble(str,first,last);

    if any( tol ~= 1 )
      xinf = dec2dble(sign,mantissa,exponent,-1);
      setround(-1)
      xinf = xinf - INTLAB_INTVAL_POWER10.sup(1,tol)';
      setround(0)
      xsup = dec2dble(sign,mantissa,exponent,+1);
      setround(1)
      xsup = xsup + INTLAB_INTVAL_POWER10.sup(1,tol)';
      setround(0)
    else
      xinf = dec2dble(sign,mantissa,exponent,-1);
      xsup = dec2dble(sign,mantissa,exponent,+1);
      setround(0)
    end
    x = infsup(xinf,xsup);

  end



function [ mant , Exp , tol ] = strmantexp2dble(str,first,last);
% String contains unsigned inputnumbers str(first(i):last(i))
% Those are split into array of mantissa digits and exponents
% output numbers .mant*10^Exp +/- 10^tol

  if isempty(first) | isempty(last) | ~isequal(size(first),size(last))
    error('Invalid input for str2intval')
  end
  n = length(first);
  L = max(last-first) + 1;                 % maximal length of inputs
  input = str2matrix(str,first,last,0)';   % inputs right bound
  index = ( input == 'e' ) | ( input == 'E' ) ;
  sumindex = sum(index,1);
  if any(sumindex>1)
    error('Invalid input for str2intval')
  end

  Exp = zeros(n,1);
  lastmant = last;
  indexexp = find(sumindex);
  if ~isempty(indexexp)
    firstexp = last(indexexp) - indexexp*L + find(index)' + 1;
    lastmant(indexexp) = firstexp - 2;
    % input(indexexp,:)  input numbers with exponent
    %   str2num(str(firstexp(i):last(indexexp(i)))) =
    %           exponents of inputnumber(indexexp(i))  for i=1:length(indexexp)

    % get signs of exponents
    nexp = length(indexexp);
    signexp = ones(nexp,1);
    v = ( str(firstexp) == '-' ) ;
    signexp(v) = -1;
    v = v | ( str(firstexp) == '+' ) ;
    firstexp(v) = firstexp(v) + 1;

    % treat exponents
    if any( last(indexexp)-firstexp < 0 )
      error('Invalid input for str2intval')
    end
    strexp = str2matrix(str,firstexp,last(indexexp),0);  % exponents right bound
    strexp(isspace(strexp)) = '0';
    if any(any( ( strexp < '0' ) | ( strexp > '9' ) ))
      error('Invalid input for str2intval')
    end

    Exp(indexexp) = signexp.*((strexp-'0')*(10.^(size(strexp,2)-1:-1:0))');
  end

  % gather mantissas
  % str2num(str(first(i):lastmant(i))) = mantissa of inputnumber(i) w/o sign
  %   for i=1:n
  index = ( input == '.' );
  sumindex = sum(index,1);
  if any(sumindex>1)
    error('An input number contains two decimal points')
  end
  % input(indexdecpt,:)  input numbers with decimal point

  indexdecpt = find(sumindex);     % treat decimal point
  if ~isempty(indexdecpt)
    decpt = zeros(n,1);
    indexpt = ( str=='.' );
    decpt(indexdecpt) = lastmant(indexdecpt) - find(indexpt);
    Exp = Exp - decpt;
    
    str(indexpt) = '';            % extract decimal points from mantissa
    index = zeros(1,n);
    index(indexdecpt) = 1;
    cumsumindex = cumsum(index);
    first = first - cumsumindex;
    first(indexdecpt) = first(indexdecpt) + 1;
    lastmant = lastmant - cumsumindex;
  end

  % strmant: mantissas w/o sign and w/o decimal point right bound
  strmant = str2matrix(str,first,lastmant,0);

  % extract tolerance '_'
  index = ( strmant == '_' );
  tol_ = sum(index,2);
  tol = ones(n,1);       % INTLAB_INTVAL_POWER10.sup(1,1) := 0 by definition
  if any(tol_)
    if any(any(diff(index,1,2)<0))     % '_' not adjacent like 3.45_45__
      error('Invalid input for str2intval')
    end
    % for indices tol_~=0: mantissas with tolerances +/- 10^(Exp+tol_)
    strmant(index) = '0';
    v = ( tol_ ~= 0 );
    tol(v) = Exp(v) + tol_(v) + 341; % INTLAB_INTVAL_POWER10.sup(1,tol) = 10^tol
  end

  strmant(isspace(strmant)) = '0';
  % strmant: mantissas w/o sign and w/o decimal point right bound,
  %   padded with 0 filling blanks

  mant = strmant - '0';
  % treat inf and NaN
  lowstrmant = lower(strmant);
  index = strmatch('fni',fliplr(lowstrmant));
  if ~isempty(index)
    mant(index,:) = repmat([zeros(1,size(mant,2)-1) 4711],length(index),1);
  end
  index = strmatch('nan',fliplr(lowstrmant));
  if ~isempty(index)
    mant(index,:) = repmat([zeros(1,size(mant,2)-1) 4711],length(index),1);
  end
  if any(any( ( ( mant < 0 ) | ( mant > 9 ) ) & ( mant~=4711 ) )) | isempty(mant)
    error('Invalid input for str2intval')
  end
  Exp = Exp+size(strmant,2)-1;


function str = extractstr(str,first,last);
% pad string with blanks except first(i):last(i)

  if any( (last-first) < 0 ) | isempty(first) | isempty(last)
    error('Invalid input to str2intval')
  end
  index = zeros(1,length(str));
  index(first) = 1;
  index(last) = -1;
  onechar = first(first==last);
  index(onechar) = 0;
  index = ( cumsum(index) == 1 );
  index(last) = 1;
  index(onechar) = 1;
  str(~index) = ' ';


function output = str2matrix(input,first,last,left)
%STR2MATRIX      Convert row string to string matrix
%
%  output = str2matrix(input,first,last,left)
%
%input    is a row string of size (1,length(input))
%           first and last are index vectors such that
%           str(i) := input(first(i):last(i))  forms a string, i=1:n
%           for n := length(first) = length(last)
%output   is a string matrix of size (n,maxlen) with maxlen=max(last-first)
%           such that the i-th row consists of str(i)
%left     if 1, str(i) is leftbound in i-th row of output, otherwise rightbound
%

  n = length(first);                    % number of inputs
  if n==1
    output = input(first:last);
  else
    len = last-first+1;                 % length of inputs excluding blanks
    L = max(len);                       % maximal length of input
    Index = zeros(1,n*L);
    if left
      i = len+L*(0:n-1);
      Index(i)=1;
      Index = ~ cumsum(reshape(Index,L,n),1);
      Index(i)=1;
    else
      i = (L-len+1)+L*(0:n-1);
      Index(i)=1;
      Index = cumsum(reshape(Index,L,n),1);
    end
    Index = Index(:);
    iota = zeros(size(input));
    iota(first) = 1;                    % vector with one's between
    iota(last) = iota(last) - 1;        % first and last
    iota = cumsum(iota);
    iota(last) = 1;
    output = blanks(n*L);
    output(logical(Index)) = input(logical(iota));
    output = reshape(output,L,n)';      % string matrix, inputs in rows
  end


function xdouble = dec2dble(s,m,e,rnd)
%DEC2DBLE     Rounded decimal to double conversion, rigorous for rnd in {-1,1}
%
%Given sign s, vector of mantissa digits m_i and exponent e,
%
%                length(m)
%                  ---        e-i+1
%  xdec  :=  s *   >   m  * 10
%                  ---  i
%                   1
%
%  call  [xinf,xsup] = dec2dble(s,m,e,rnd)
%
%  rnd = -1 :  rounding downwards, xdouble <= xdec
%         1                                   xdec <= xdouble
%         0 :  approximate rounding to nearest
%
%For m being a matrix, mantissas are rows of m,
%  for example
%    s = [  1 ; ...
%          -1 ]
%    m = [  3 1 4 1 5 9 ; ...
%           2 7 1 8 2 8 ] ;
%    e = [  0 ; ...
%           0 ]
%
% elements of e in [-309,309] ,   elements of m in [0,9]
%

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');

  offset = 9*(e+341)*ones(1,size(m,2));
  mantissa = ones(size(m,1),1)*(9*(size(m,2):-1:1));

  infindex = ( m(:,end)==4711 );
  if ~isempty(infindex)
    m(infindex,end) = 0;
  end
  pos = ( s==1 );
  mpos= fliplr(m(pos,:));
  % treat inf
  posindex = (mpos~=0).*(offset(pos,:) - mantissa(pos,:) + mpos);
  posindex = max(posindex,2) - (mpos==0);
  posindex(posindex>5841) = 5842;
  neg = ~pos;
  mneg = fliplr(m(neg,:));
  negindex = (mneg~=0).*(offset(neg,:) - mantissa(neg,:) + mneg);
  negindex = max(negindex,2) - (mneg==0);
  negindex(negindex>5841) = 5842;

  if rnd==-1
    xdouble = zeros(size(m,1),1);
    setround(-1)
    xdouble(pos) = sum( INTLAB_INTVAL_POWER10.inf(posindex) , 2 );
    setround(1)
    xdouble(neg) = - sum( INTLAB_INTVAL_POWER10.sup(negindex) , 2 );
  elseif rnd==1
    xdouble = zeros(size(m,1),1);
    setround(1)
    xdouble(pos) = sum( INTLAB_INTVAL_POWER10.sup(posindex) , 2 );
    setround(-1)
    xdouble(neg) = - sum( INTLAB_INTVAL_POWER10.inf(negindex) , 2 );
  elseif rnd==0
    xinf = zeros(size(m,1),1);
    xsup = xinf;
    xinf(pos) = sum( INTLAB_INTVAL_POWER10.inf(posindex) , 2 );
    xinf(neg) = - sum( INTLAB_INTVAL_POWER10.sup(negindex) , 2 );
    xsup(pos) = sum( INTLAB_INTVAL_POWER10.sup(posindex) , 2 );
    xsup(neg) = - sum( INTLAB_INTVAL_POWER10.inf(negindex) , 2 );
    xdouble = xinf + 0.5*(xsup-xinf);
  else
    error('invalid call of dec2double')
  end
  if ~isempty(infindex)
    xdouble(infindex) = NaN;
  end
