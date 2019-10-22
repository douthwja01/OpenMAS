function outstr = disp_(c,name,restricted)
%DISP_        Display interval with uncertainty (rigorous)
%
%The call
%
%   disp_(c)        diplays interval c using underscore for uncertainty digits
%   str = disp_(c)  puts output into string str
%
%using Matlab formats "short (e)" and "long (e)".  Correct interval is always
%  last displayed figure +/- 1 .
%Output in string str columnwise; can be used for input by intval(str).
%  For safety, output string is appended with '_' to indicate interval.
%
%For example,  infsup(3.14159,3.14160)  produces in "format short"
%  intval ans =
%     3.1416
%which is [3.1415,3.1417]. Wide intervals like infsup(-12,-7) produce
%  intval ans =
%   -1_.____
%which is in fact [0,20], and very wide infsup(-3,8) produce
%  intval ans =
%    0_.____
%which is [-10,10]. Complex intervals are displayed the same way, for example
%midrad(4711-.1i,1e-10) produces in "format long"
%  intval ans =
%   1.0e+003 *
%    4.711000000000__ -  0.000100000000__i
%and inf "format long e"
%  intval ans =
%   4.711000000000___e+003 - 1.00000000_______e-001i
%
%If intervals are too wide, the display is automatically changed
%to infimum/supremum for real intervals and to midpoint/radius for
%complex intervals.
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
% modified 06/23/98     S.M. Rump  multi-dimensional arrays, improved display,
%                                  input near realmax, sparse arrays
% modified 08/29/00     S.M. Rump  output value zero for output string with _
%                                  output for intval(-inf)
%                                  special output in outstr.exp and .str added
%                                  rounding unchanged after use
%
% modified 10/03/02     S.M. Rump  improved for Matlab 6.1f and sparse output
% modified 01/18/04     S.M. Rump  argument nargout_ in intval2str instead of ambiguous nargout,
%                                    rounding unaltered after return, output for IBM corrected
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    infinity corrected
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%                                    linear indices for huge matrices
% modified 01/03/07     S.M. Rump  Automatic change of output if IVs too wide
% modified 05/22/07     S.M. Rump  test for mean([]) 
% modified 08/25/07     S.M. Rump  huge indices for sparse matrices
% modified 02/17/08     S.M. Rump  Matlab bug fix
% modified 10/14/08     S.M. Rump  huge indices for sparse matrices, factor removed
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/04/12     S.M. Rump  take care of sparse input
%

  INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
  if isequal(INTLAB_INTVAL_DISPLAY,'Display_')   % otherwise forced to use disp_
    wng = warning;
    warning off
%VVVV nnzc = nonzeros(c);
    [dummy1,dummy2,nnzc] = find(c);
%AAAA Matlab 2013b bug fix: " a=intval(0) " causes crash and core dump
%VVVV M = relerr(nnzc(:));
      s.type = '()';
      s.subs = {':'};
      M = relerr(subsref(nnzc,s));
%AAAA Matlab V5.2 bug fix
    M = mean(M(isfinite(M(:))));
    warning(wng)
    if ( ~isempty(c) ) & ( M>.02 )   % intervals inaccurate, change to other output
      if nargout
        evalstr = 'outstr = '
      else
        evalstr = '';
      end
      if c.complex
        evalstr = [ evalstr 'midrad(c' ];
      else
        evalstr = [ evalstr 'infsup(c' ];
      end
      if nargin==1
        eval([ evalstr ','''',1);' ]);
      else
        evalstr = [ evalstr ',name' ];
        if nargin>=3                      % this must be true
          evalstr = [ evalstr ',restricted' ];
        end
        eval([ evalstr ');' ]);
      end
      return
    end
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==1
    name = inputname(1);
    restricted = 0;
  elseif nargin==2
    error('invalid call of disp_')
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
      disp_(cs,[name '(:,:' str ')'],0);
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

  commonexp = 0;
  if ~isequal(format(end),'E')
    % calculate exponent range and common factor
    if c.complex
      % cmid(:)+crad(:) surprisingly slow for sparse c
      if isequal(crad,0)
        cre = abs(real(cmid));
        cim = abs(imag(cmid));
      else
        cre = abs(real(cmid)) + crad;
        cim = abs(imag(cmid)) + crad;
      end
      cre = nonzeros(cre);
      cim = nonzeros(cim);
      cre(isinf(cre)) = 0;                      % fast for sparse matrices
      cim(isinf(cim)) = 0;
      cremax = max(cre(:));
      cimmax = max(cim(:));
      cmax = max(cremax,cimmax);
    else
      if issparse(crad) & ~any(any(crad~=0))
        cmax = abs(cmid);
      else
        cmax = abs(cmid) + crad;
      end
      cmax = nonzeros(cmax);
      cmax(isinf(cmax)) = 0;                    % fast for sparse matrices
      cmax = max(cmax(:));
    end
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
  else
    switch format
      case 'short',  len = 10; prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
  end
  INTLAB_DISPLAY_WIDTH = getappdata(0,'INTLAB_DISPLAY_WIDTH');    % minimum 110, so columns>=1
  columns = floor(INTLAB_DISPLAY_WIDTH/len);

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
      outstr.str = char(zeros(m,2*len+3));
      for i=1:m
        s = '';
        chsign = '+';
        if imag(cmid(i)) < 0
          chsign = '-';
        end
        if isequal(crad,0)
          cradij = 0;
        else
          cradij = crad(i);
        end
        cmidi = full(cmid(i));
        cradij = full(cradij);
        strimag = intval2str(abs(imag(cmidi)),cradij, ...
                             commonexp,len,prec,expon,0);
        str = intval2str(real(cmidi),cradij,commonexp,len,prec,expon,0);
        s = [ s str ' ' chsign strimag(2:end) 'i ' ];
        outstr.str(i,:) = s;
      end
    else                   % real intervals
      outstr.str = char(zeros(m,len));
      for i=1:m
        s = intval2str(full(cmid(i)),full(crad(i)),commonexp,len,prec,expon,0);
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
        if isequal(crad,0)
          cradi = 0;
        else
          cradi = crad(i);
        end
        cmidi = full(cmid(i));
        cradi = full(cradi);
        outstr = [ outstr ...
                   intval2str(real(cmidi),cradi,commonexp,len,prec,expon,1) ];
        strim = intval2str(imag(cmidi),cradi,commonexp,len,prec,expon,1);
        strim(isspace(strim)) = '';      % eliminate spaces
        if strim(1)~='-'
          strim = [ '+' strim ];
        end
        outstr = [ outstr strim 'i' ];
      end
    else
      for i=1:prod(size(cmid))
        outstr = [ outstr ...
                   intval2str(full(cmid(i)),full(crad(i)),commonexp,len,prec,expon,1) ];
      end
    end
    setround(rndold)
    return
  end

  if issparse(c)
    if isequal(crad,0)                              % mid or rad maybe zero
      [I,J] = find(spones(cmid));    
    else
      [I,J] = find(spones(cmid)+spones(crad));  
    end
    if length(I)==0
      dummy = mid(c);
      disp(dummy)       % make sure empty matrices are displayed
      setround(rndold)
      return
    end
    if c.complex
      for i=1:length(I)
        if isequal(crad,0)
          cradi = 0;
        else
          cradi = crad(I(i),J(i));
        end
        cmidij = full(cmid(I(i),J(i)));
        cradi = full(cradi);
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        str = [ str intval2str(real(cmidij),cradi,commonexp,len,prec,expon,0) ];
        strim = intval2str(abs(imag(cmidij)),cradi,commonexp,len,prec,expon,0);
        chsign = '+';
        if imag(cmidij) < 0
          chsign = '-';
        end
        str = [ str ' ' chsign strim(3:end) 'i ' ];
        disp( str )
      end
    else
      for i=1:length(I)
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        disp([str ...
        intval2str(full(cmid(I(i),J(i))),full(crad(I(i),J(i))),commonexp,len,prec,expon,0) ])
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
          chsign = '+';
          if imag(cmid(i,j)) < 0
            chsign = '-';
          end
          if isequal(crad,0)
            cradij = 0;
          else
            cradij = crad(i,j);
          end
          cmidij = full(cmid(i,j));
          cradij = full(cradij);
          strimag = intval2str(abs(imag(cmidij)),cradij, ...
                               commonexp,len,prec,expon,0);
          str = intval2str(real(cmidij),cradij,commonexp,len,prec,expon,0);
          s = [ s str ' ' chsign strimag(2:end) 'i ' ];
        end
        disp(s)
      end
    else                   % real intervals
      for i=1:m
        s = '';
        for j = j1:j2
          str = intval2str(full(cmid(i,j)),full(crad(i,j)),commonexp,len,prec,expon,0);
          s = [ s str ];
        end
        disp(s)
      end
    end
    if loose & ~restricted, disp(' '), end
  end
    
  setround(rndold)


function str = intval2str(cmid,crad,commonexp,len,prec,expon,nargout_)
%INTVAL2STR   internal routine for output routine display
%

%Rigorous conversion of interval  cmid+/-crad  with directed rounding
%  and common exponent; only for one real input.
%Output format is %len.precX with  X=e/f  for  expon=1/0, respectively.
%In case expon=0 (i.e. f-format), 0 <= c <= 1000 required (format short or long)
%  commonexp~=0 only in case expon=0
%For output format: 2 <= prec <= 15 with len large enough to carry all digits
%
%String "str" is constructed according to format such that content of
%  string "str" contains  cmid+/-crad
%
%If nargout_, extra '_' added to indicate interval
%

  cmid = full(cmid);
  crad = full(crad);
  % exception handling
  if isnan(cmid) | isnan(crad)
    str = [ blanks(len-3) 'NaN' ];
    return
  end
  if isinf(cmid) | isinf(crad)
    if isinf(crad)
      str = [ blanks(len-6) '+/-Inf' ];
    else
      if cmid<0
        str = [ blanks(len-4) '-Inf' ];
      else
        str = [ blanks(len-3) 'Inf' ];
      end
    end
    return
  end

  cmidrnd = cmid;
  format = [ '%' sprintf('%d',len) '.' sprintf('%d',prec) 'e' ];
  if ~expon
    format(end) = 'f';
  end

  if ( cmid == 0 ) & ( crad == 0 )    % interval 0
    str = sprintf(format,0);
    if nargout_
      if expon
        str = [ str(1:end-5) '_' str(end-4:end) ];
      else
        str = [ str '_' ];
      end
    end
    return
  end

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');

  cmidrnd = abs(cmidrnd);
  signcmid = sign(cmid);

  % compute internal display format, "pprec" figures after decimal point
  factor = INTLAB_INTVAL_POWER10.sup(1,commonexp+341);   % 10^commonexp
  ws = warning;
  warning off
  e1 = floor(log10(cmidrnd));
  e2 = floor(log10(crad));
  warning(ws)

  precmax = prec;
  if ~expon & ~commonexp
    precmax = prec+e1;
  end
  pprec = min(e1-e2-1,precmax);
  if ~expon
    pprec = pprec+1;     % improved output for midrad(1-1e-8,1e-14)
  end


  % 3 stages of printing precision (digits after decimal point) and output:
  %
  % stage 1:  e1 > e2  radius small compared to midpoint, cmidrnd = abs(cmid)
  %                    printing precision pprec >= 0
  % stage 2:  e1 = e2  radius of same size as midpoint, pprec < 0
  %     2.1:  crad <= abs(cmid)  sign determined, cmidrnd = 10^(e1+1)
  %                              wide = 0
  %     2.2:  crad > abs(cmid)   zero interval, cmidrnd = 1 * 10^(e2+1)
  %                              wide = 1, leading 1 is padded with 0
  %           e1 < e2  same as 2.2
  %
  wide = 0;
  if pprec<1
    pprec = 0;
    if crad <= abs(cmid)
      cmidrnd = INTLAB_INTVAL_POWER10.sup(1,e1+1+341);   % 10^(e1+1)
    else
      wide = 1;
      cmidrnd = INTLAB_INTVAL_POWER10.sup(1,e2+1+341);   % 10^(e2+1)
    end
  end
  if isinf(cmidrnd)
    cmidrnd = sign(cmidrnd) * realmax;
  end

  while 1
%AAAA bug fix for Watcom underflow printf
    if (cmidrnd>0) & (cmidrnd<realmin)
      formatstr = '%22.15e';
%VVVV end bug fix
    else
      formatstr = [ '%' sprintf('%d',pprec+7) '.' sprintf('%d',pprec) 'e' ];
    end
    str = sprintf(formatstr,cmidrnd/factor);
    % make sure that exponent has three digits
    if ( lower(str(end-3))=='e' )
      str = [ str(2:end-2) '0' str(end-1:end) ];
    end

    % convert str "x.yyyyy_" with decimal point after x and "pprec" y's and
    %   uncertainty in the last displayed figure back to intval [xinf,xsup]
    %   with correct rounding
    indexdp = find( str=='.' );
    if isempty(indexdp)
      str = str(2:end);        % one output digit, no decimal point
    else
      str( str=='.' ) = '';
    end
    if wide
      str(1) = '0';
    end

    Exp = ( 2*(str(end-3)=='+') - 1 ) * (str(end-2:end)-'0')*[100;10;1];
    indexend = pprec + 1;
    Exp_ = commonexp+Exp-indexend+1;   % exponent of uncertainty

    m = str(indexend:-1:1) - '0';

    % convert back to double
    offset = 9*(Exp+commonexp+341);
    mant = 9*(indexend:-1:1);
    index = (m~=0).*(offset - mant + m) + (m==0);
    
    if signcmid==1
      setround(-1)
      xinf = sum( INTLAB_INTVAL_POWER10.inf(index) ) - ...
                  INTLAB_INTVAL_POWER10.sup(1,Exp_+341) ;
      setround(1)
      xsup = sum( INTLAB_INTVAL_POWER10.sup(index) ) + ...
                  INTLAB_INTVAL_POWER10.sup(1,Exp_+341) ;
    else
      setround(-1)
      xsup = - ( sum( INTLAB_INTVAL_POWER10.inf(index) ) - ...
                      INTLAB_INTVAL_POWER10.sup(1,Exp_+341) );
      setround(1)
      xinf = - ( sum( INTLAB_INTVAL_POWER10.sup(index) ) + ...
                      INTLAB_INTVAL_POWER10.sup(1,Exp_+341) );
    end
    cont = ( xinf+crad <= cmid ) & ( cmid+crad <= xsup );
    setround(0)
    
    if cont             % output not yet valid
      % current string correct: "x.yyyyy" with decimal point after x
      % and "pprec" y's and uncertainty in the last displayed figure.
      % Convert that string to desired format %len.precX with
      %   X = e/f according to expon = 1/0, respectively.
      if expon     % e-format
        str = [ blanks(len-prec-7) str(1) '.' str(2:pprec+1) ...
                char('_'+zeros(1,prec-pprec)) str(end-4:end) ];
        if signcmid == -1
          str(len-prec-7) = '-';
        end
      else                      % f-format
        digdp = pprec - Exp;    % digits after decimal point
        k = pprec+1 - max(digdp-prec,0);  % digits to be used
        mantissa = [ char('0'+zeros(1,min(-Exp,prec+1))) ...
                     str(1:k) char('_'+zeros(1,prec-digdp)) ];
        if signcmid == -1
          mantissa = [ '-' mantissa ];
        end
        str = [ blanks(len-length(mantissa)-1) mantissa ];
        str = [ str(1:len-prec-1) '.' str(len-prec:end) ];
      end
      break
    end

    if pprec == 0
      cmidrnd = INTLAB_INTVAL_POWER10.sup(1,Exp_+1+341);
    else
      pprec = pprec - 1;
    end

  end

  if nargout_
    if expon
      str = [ str(1:end-5) '_' str(end-4:end) ];
    else
      str = [ str '_' ];
    end
  end
