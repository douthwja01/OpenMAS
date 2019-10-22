function outstr = midrad(c,name,restricted)
%MIDRAD       Display of intervals by midpoint/radius (rigorous)
%
%The call
%
%  midrad(c)        displays interval c in < mid , rad > representation
%  str = midrad(c)  puts output into string str
%
%Output in string str columnwise; can be used for input by intval(str).
%

%for internal use:
%  name                name of output variable
%  restricted == 1     no header, no extra lines output
%
%Call only with 1 or 3 input arguments
%
%Special call:  outstr = midrad(x,[],[])
%  for column vector x output in outstr.exp and outstr.str
%

% written  11/30/98     S.M. Rump
% modified 06/10/98     S.M. Rump  multi-dimensional arrays, exceptions
% modified 06/23/99     S.M. Rump  dble2str -> dble2str_rnd in \private,
%                                  sparse arrays
% modified 08/29/00     S.M. Rump  special output in outstr.exp and .str added
%                                  rounding unchanged after use
% modified 10/03/02     S.M. Rump  improved for Matlab 6.1f and sparse output
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%                                    linear indices for huge matrices
% modified 08/25/07     S.M. Rump  huge indices for sparse matrices
% modified 10/14/08     S.M. Rump  huge indices for sparse matrices
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/04/12     S.M. Rump  take care of sparse input
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin<=2
    name = inputname(1);
    restricted = 0;
  end

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  siz = size(c);
  if ~restricted & ~nargout & ( length(siz)<3 )
    if loose, disp(' '); end
    disp([ 'intval ' name ' = ' ])
    if loose, disp(' '); end
  end

  if isempty(c)          % empty interval
    if c.complex
      c.mid
    else
      c.inf
    end
    if loose & ~restricted & ~nargout
      disp(' ')
    end
    setround(rndold)
    return
  end

  if length(siz)>2 & ~nargout
    siz = siz(3:end);
    len = length(siz);
    prodsiz = cumprod([1 siz(1:len-1)]);

    for k=0:prod(siz)-1
      y = mod(floor(k./prodsiz),siz);
      str = sprintf(',%d',y+1);
%VVVV eval(['cs = c(:,:' str ');']);
      s.type = '()';
      eval(['s.subs = {'':'','':''' str '};']);
      cs = subsref(c,s);
%AAAA Matlab V5.2 bug fix
      midrad(cs,[name '(:,:' str ')'],0);
    end
    setround(rndold)
    return
  end

  cmid = mid(c);
  crad = rad(c);
  [m n] = size(cmid);

  format = get(0,'format');
  if ~isequal(format,'short') & ~isequal(format,'long') & ...
     ~isequal(format,'shortE') & ~isequal(format,'longE')
    format = 'short';
  end

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');
  commonexp = 0;
  if ~isequal(format(end),'E')
    % calculate exponent range and common factor
    if c.complex
      cre = abs(real(cmid));
      cim = abs(imag(cmid));
      cre = nonzeros(cre);
      cim = nonzeros(cim);
      if isempty(cre)                   % careful: max(1,[]) = []
        cre = 0;
      end
      if isempty(cim)
        cim = 0;
      end
      cre(isinf(cre)) = 0;
      cim(isinf(cim)) = 0;
      cremax = max(cre(:));
      cimmax = max(cim(:));
      cmax = max(cremax,cimmax);
    else
      cmid_ = abs(cmid);
      cmid_ = nonzeros(cmid_);
      if isempty(cmid_)                 % careful: max(1,[]) = []
        cmid_ = 0;
      end
      cmid_(isinf(cmid_)) = 0;
      cmax = max(cmid_);
      cmax = max(cmax(:));
    end
    crad_ = abs(crad);
    crad_ = nonzeros(crad_);
      if isempty(crad_)                 % careful: max(1,[]) = []
        crad_ = 0;
      end
    crad_(isinf(crad_)) = 0;
    cmax = max( cmax , max(crad_(:)) );
    % care for sparse matrices
    cmax = full(cmax);
    if isempty(cmax)
      cmax = 0;
    end

    emax = floor(log10( cmax + (cmax==0) ));
    if c.complex
      if emax >= 2
        commonexp = emax;
      end
      if emax <= -3
        commonexp = emax+1;
      end
    else
      if isequal(format,'short') & emax >= 3
        commonexp = emax;
      elseif isequal(format,'long') & emax >= 2
        commonexp = emax;
      end
      if emax <= -4
        commonexp = emax+1;
      end
    end
  end
  factor = INTLAB_INTVAL_POWER10.sup(1,commonexp+341);    % 10^commonexp

  if commonexp~=0                      % common factor
    if nargout                         % output to str: be sure with exponent
      if ~( isempty(name) & isempty(restricted) )   % if not special output
        commonexp = 0;
        if ~isequal(format(end),'E')
          format = [ format 'E' ];
        end
      end
    else
      strexp = sprintf('%04d',commonexp);
      if commonexp >= 0
        strexp(1) = '+';
      end
      disp([ '  1.0e' strexp ' *' ])
      if loose, disp(' '); end
    end
  end

  if c.complex
    switch format
      case 'short',  len = 9;  prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
    len1 = 3*len+5;    % length of one element in current format
  else
    switch format
      case 'short',  len = 10; prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
    len1 = 2*len+2;    % length of one element in current format
  end
  INTLAB_DISPLAY_WIDTH = getappdata(0,'INTLAB_DISPLAY_WIDTH');    % minimum 110, so columns>=1
  columns = floor((INTLAB_DISPLAY_WIDTH+1)/(len1+1));

  formatstr = [ '%' int2str(len) '.' int2str(prec) 'f' ];
  if expon
    formatstr(end) = 'e';
  end

  if nargout & isempty(name) & isempty(restricted)
    % special case: input c column vector, output in outstr.exp and outstr.str

    if commonexp~=0                      % common factor
      strexp = sprintf('%04d',commonexp);
      if commonexp >= 0
        strexp(1) = '+';
      end
      outstr.exp = [ '  1.0e' strexp ' *' ];
    else
      outstr.exp = '';
    end

    if c.complex           % complex intervals
      outstr.str = char(zeros(m,3*len+5));
      for i=1:m
        s = '';
        cmidi = full(cmid(i));
        strrealmid = sprintf(formatstr,real(cmidi)/factor);
        signimagmid = sign(imag(cmidi));
        strimagmid = sprintf(formatstr,abs(imag(cmidi))/factor);
        strimagmid = strimagmid(2:end);     % omit leading space
        [cmidrealinf,cmidrealsup] = str2intval(strrealmid,commonexp);
        if signimagmid==1
          [cmidimaginf,cmidimagsup] = str2intval(strimagmid,commonexp);
        else
          [cmidimagsup,cmidimaginf] = str2intval(strimagmid,commonexp);
          cmidimagsup = - cmidimagsup ;
          cmidimaginf = - cmidimaginf ;
        end
        setround(1)
        re = max( real(cmidi)-cmidrealinf,cmidrealsup-real(cmidi) );
        im = max( imag(cmidi)-cmidimaginf,cmidimagsup-imag(cmidi) );
        if isequal(crad,0)
          cradij = abs( re + sqrt(-1)*im );
        else
          cradij = crad(i) + abs( re + sqrt(-1)*im );
        end
        cradij = full(cradij);
        chsign = '+';
        if signimagmid==-1
          chsign = '-';
        end
        strrad = dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1);
        s = [ s '<' strrealmid ' ' chsign strimagmid 'i,' ...
              strrad '> ' ];
        outstr.str(i,:) = s;
      end
    else                   % real intervals
      outstr.str = char(zeros(m,2*len+3));
      for i=1:m
        s = '';
        cmidi = full(cmid(i));
        strmid = sprintf(formatstr,cmidi/factor);
        [cmidinf,cmidsup] = str2intval(strmid,commonexp);
        setround(1)
        cradij = full( crad(i) + max(cmidi-cmidinf,cmidsup-cmidi) );
        s = [ s '<' strmid ',' ...
              dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1) '> ' ];
        outstr.str(i,:) = s;
      end
    end
    setround(rndold)
    return
  end

  if nargout
    outstr = [];
    if c.complex
      for i=1:prod(size(cmid))
        cmidi = full(cmid(i));
        strrealmid = sprintf(formatstr,real(cmidi));
        strimagmid = sprintf(formatstr,imag(cmidi));
        strimagmid(isspace(strimagmid)) = '';
        if strimagmid(1)~='-'
          strimagmid = [ '+' strimagmid ];
        end
        [cmidrealinf,cmidrealsup] = str2intval(strrealmid,0);
        [cmidimaginf,cmidimagsup] = str2intval(strimagmid,0);
        setround(1)
        re = max( real(cmidi)-cmidrealinf,cmidrealsup-real(cmidi) );
        im = max( imag(cmidi)-cmidimaginf,cmidimagsup-imag(cmidi) );
        if isequal(crad,0)
          cradi = abs( re + sqrt(-1)*im );
        else
          cradi = crad(i) + abs( re + sqrt(-1)*im );
        end
        cradi = full(cradi);
        strrad = dble2str_rnd(cradi,commonexp,len-1,prec,expon,+1);
        outstr = [ outstr '<' strrealmid strimagmid 'i,' strrad '> ' ];
      end
    else
      for i=1:prod(size(cmid))
        cmidi = full(cmid(i));
        strmid = sprintf(formatstr,cmidi);
        [cmidinf,cmidsup] = str2intval(strmid,0);
        setround(1)
        cradi = full( crad(i) + max(cmidi-cmidinf,cmidsup-cmidi) );
        outstr = [ outstr '<' strmid ',' ...
              dble2str_rnd(cradi,commonexp,len-1,prec,expon,+1) '> ' ];
      end
    end
    setround(rndold)
    return
  end

  if issparse(c)
    if isequal(crad,0)                            % mid or rad maybe zero
      [I,J] = find(spones(cmid));   
    else
      [I,J] = find(spones(cmid)+spones(crad));    
    end
    if length(I)==0
      mid(c)
      setround(rndold)
      return
    end
    if c.complex
      for i=1:length(I)
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        cmidij = full(cmid(I(i),J(i)));
        strrealmid = sprintf(formatstr,real(cmidij)/factor);
        signimagmid = sign(imag(cmidij));
        strimagmid = sprintf(formatstr,abs(imag(cmidij))/factor);
        strimagmid = strimagmid(2:end);     % omit leading space
        [cmidrealinf,cmidrealsup] = str2intval(strrealmid,commonexp);
        if signimagmid==1
          [cmidimaginf,cmidimagsup] = str2intval(strimagmid,commonexp);
        else
          [cmidimagsup,cmidimaginf] = str2intval(strimagmid,commonexp);
          cmidimagsup = - cmidimagsup ;
          cmidimaginf = - cmidimaginf ;
        end
        setround(1)
        re = max( real(cmidij)-cmidrealinf,cmidrealsup-real(cmidij) );
        im = max( imag(cmidij)-cmidimaginf,cmidimagsup-imag(cmidij) );
        if isequal(crad,0)
          cradij = abs( re + sqrt(-1)*im );
        else
          cradij = crad(I(i),J(i)) + abs( re + sqrt(-1)*im );
        end
        cradij = full(cradij);
        chsign = '+';
        if signimagmid==-1
          chsign = '-';
        end
        strrad = dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1);
        str = [ str '<' strrealmid ' ' chsign strimagmid 'i,' ...
                strrad '> ' ];
        disp( str )
      end
    else
      for i=1:length(I)
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        cmidij = full(cmid(I(i),J(i)));
        strmid = sprintf(formatstr,cmidij/factor);
        [cmidinf,cmidsup] = str2intval(strmid,commonexp);
        setround(1)
        cradij = full(crad(I(i),J(i)) + max(cmidij-cmidinf,cmidsup-cmidij));
        str = [ str '<' strmid ',' ...
                dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1) '> ' ];
        disp( str )
      end
    end
    setround(rndold)
    return
  end

  for jj=1:ceil(n/columns)
    j1 = (jj-1)*columns+1;
    if jj*columns<n
      j2 = jj*columns;
    else
      j2 = n;
    end
    if n>columns
      if j1~=j2
        disp(['  Columns ' sprintf('%d',j1) ' through ' sprintf('%d',j2)]);
      else
        disp(['  Column ' sprintf('%d',j1)]);
      end
      if loose, disp(' '); end
    end
    if c.complex           % complex intervals
      for i=1:m
        s = '';
        for j = j1:j2
          cmidij = full(cmid(i,j));
          strrealmid = sprintf(formatstr,real(cmidij)/factor);
          signimagmid = sign(imag(cmidij));
          strimagmid = sprintf(formatstr,abs(imag(cmidij))/factor);
          strimagmid = strimagmid(2:end);     % omit leading space
          [cmidrealinf,cmidrealsup] = str2intval(strrealmid,commonexp);
          if signimagmid==1
            [cmidimaginf,cmidimagsup] = str2intval(strimagmid,commonexp);
          else
            [cmidimagsup,cmidimaginf] = str2intval(strimagmid,commonexp);
            cmidimagsup = - cmidimagsup ;
            cmidimaginf = - cmidimaginf ;
          end
          setround(1)
          re = max( real(cmidij)-cmidrealinf,cmidrealsup-real(cmidij) );
          im = max( imag(cmidij)-cmidimaginf,cmidimagsup-imag(cmidij) );
          if isequal(crad,0)
            cradij = abs( re + sqrt(-1)*im );
          else
            cradij = crad(i,j) + abs( re + sqrt(-1)*im );
          end
          cradij = full(cradij);
          chsign = '+';
          if signimagmid==-1
            chsign = '-';
          end
          strrad = dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1);
          s = [ s '<' strrealmid ' ' chsign strimagmid 'i,' ...
                strrad '> ' ];
        end
        disp(s)
      end
    else                   % real intervals
      for i=1:m
        s = '';
        for j = j1:j2
          cmidij = full(cmid(i,j));
          strmid = sprintf(formatstr,cmidij/factor);
          [cmidinf,cmidsup] = str2intval(strmid,commonexp);
          setround(1)
          cradij = full( crad(i,j) + max(cmidij-cmidinf,cmidsup-cmidij) );
          s = [ s '<' strmid ',' ...
                dble2str_rnd(cradij,commonexp,len-1,prec,expon,+1) '> ' ];
        end
        disp(s)
      end
    end
    if loose & ~restricted, disp(' '); end
  end

  setround(rndold)


function [xinf,xsup] = str2intval(str,commonexp)
%STR2INTVAL   internal routine for output routines midrad
%

%Rigorous conversion of string to lower/upper bound with directed rounding
%  and common exponent; only for one real input.
%Input string in e- or f-format, possibly with exponent
%Input number is assumed to be correct.
%

  % exception handling
  if any( str=='a' ) | any( str=='f' )    % NaN or Inf
    xinf = NaN;
    xsup = NaN;
    return
  end

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');
  [m i] = max(str~=' ');
  str(1:i-1) = '';             % omit leading blanks

  sign = 1;                    % calculate sign
  if str(1)=='+'
    str = str(2:end);
  elseif str(1)=='-'
    sign = -1;
    str = str(2:end);
  end

  indexdp = find(str=='.');    % index of decimal point
  str(indexdp) = '';

  indexexp = find(str=='e');   % check exponent
  if ~isempty(indexexp)
    Exp = ( 2*(str(indexexp+1)=='+') - 1 ) * ...
          (str(indexexp+2:end)-'0')*(10.^((length(str)-indexexp-2):-1:0)');
    indexend = indexexp-1;
  else
    Exp = 0;
    indexend = length(str);
  end

  m = str(indexend:-1:1) - '0';

  % convert back to double
  offset = 9*(indexdp-2+Exp+commonexp+341);
  mant = 9*(indexend:-1:1);
  index = (m~=0).*(offset - mant + m) + (m==0);
  if sign==1
    setround(-1)
    xinf = sum( INTLAB_INTVAL_POWER10.inf(index) );
    setround(1)
    xsup = sum( INTLAB_INTVAL_POWER10.sup(index) );
  else
    setround(1)
    xinf = - sum( INTLAB_INTVAL_POWER10.sup(index) );
    setround(-1)
    xsup = - sum( INTLAB_INTVAL_POWER10.inf(index) );
  end
  setround(0)
