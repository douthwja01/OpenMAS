function outstr = infsup(c,name,restricted)
%INFSUP       Display of intervals by infimum/supremum (rigorous)
%
%The call
%
%  infsup(c)        displays interval c in [ inf , sup ] representation
%  str = infsup(c)  puts output into string str
%
%Output in string str columnwise; can be used for input by intval(str).
%

%for internal use:
%  name                name of output variable
%  restricted == 1     no header, no extra lines output
%
%Call only with 1 or 3 input arguments
%
%Special call:  outstr = infsup(x,[],[])
%  for column vector x output in outstr.exp and outstr.str
%

% written  11/30/98     S.M. Rump
% modified 06/10/98     S.M. Rump  multi-dimensional arrays, exceptions
% modified 06/23/99     S.M. Rump  dble2str -> dble2str_rnd in \private,
%                                  sparse arrays
% modified 08/29/00     S.M. Rump  special output in outstr.exp and .str added
%                                  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 08/25/07     S.M. Rump  huge indices for sparse matrices
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/14/08     S.M. Rump  huge indices for sparse matrices
% modified 06/15/11     S.M. Rump  allow inner inclusions
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

  % a=sparse([],[],[]), [i,j,s]=find(a), b=a(i), c=min(a,b), islogical(c)
  % gives 0 in version 6.0, but 1 in version 6.5 and 6.5.1
  if islogical(c.inf) 
    c.inf = double(c.inf);
  end
  if islogical(c.sup) 
    c.sup = double(c.sup);
  end
  
  if nargin==2
    if isequal(name,'inner')
      InnerInclusion = -1;
    else
      error('invalid call of infsup')
    end
  else
    InnerInclusion = 1;
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
    if rndold
      setround(rndold)
    end
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
      infsup(cs,[name '(:,:' str ')'],0);
    end
    if rndold
      setround(rndold)
    end
    return
  end

  cinf = inf(c);
  csup = sup(c);
  [m n] = size(cinf);

  format = get(0,'format');
  if ~isequal(format,'short') & ~isequal(format,'long') & ...
     ~isequal(format,'shortE') & ~isequal(format,'longE')
    format = 'short';
  end

  commonexp = 0;
  if ~isequal(format(end),'E')
    % calculate exponent range and common factor
    if c.complex
      cinf_ = abs(nonzeros(cinf));
      csup_ = abs(nonzeros(csup));
      if isempty(cinf_)                      % careful: max(1,[]) = [] 
        cinf_ = 0;
      end
      if isempty(csup_)
        csup_ = 0;
      end
      cinf_(isinf(cinf_)) = 0;
      csup_(isinf(csup_)) = 0;
      cinf_max = max(cinf_(:));
      csup_max = max(csup_(:));
      cmax = max( [ cinf_max csup_max ] );   % careful: max(1,[]) = []
    else
      cinf_ = abs(nonzeros(cinf));
      csup_ = abs(nonzeros(csup));
      if isempty(cinf)                       % careful: max(1,[]) = []
        cinf = 0;
      end
      if isempty(csup)
        csup = 0;
      end
      cinf_(isinf(cinf_)) = 0;
      csup_(isinf(csup_)) = 0;
      cinf_max = max(cinf_(:));
      csup_max = max(csup_(:));
      cmax = max( [ cinf_max csup_max ] );   % careful: max(1,[]) = []
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
      if loose, disp(' '), end
    end
  end

  if c.complex
    switch format
      case 'short',  len = 9;  prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
    len1 = 4*len+7;    % length of one element in current format
  else
    switch format
      case 'short',  len = 10; prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
    len1 = 2*len+3;    % length of one element in current format
  end
  INTLAB_DISPLAY_WIDTH = getappdata(0,'INTLAB_DISPLAY_WIDTH');    % minimum 110, so columns>=1
  columns = floor((INTLAB_DISPLAY_WIDTH+1)/(len1+1));

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
      outstr.str = char(zeros(m,4*len+8));
      for i=1:m
        s = '';
        s = [ s '[' ];
        for ii=[-1 1]
          if ii==-1
            x = full(cinf(i));
          else
            x = full(csup(i));
          end
          strre = dble2str_rnd(real(x),commonexp,len,prec,expon,ii);
          if imag(x)==0
            chsign = ' ';
            strim = blanks(len);
          else
            signim = sign(imag(x));
            chsign = '+';
            if signim==-1
              chsign = '-';
            end
            strim = [ dble2str_rnd(abs(imag(x)), ...
                      commonexp,len-1,prec,expon,ii*signim) 'i' ];
          end
          s = [ s strre ' ' chsign strim ];
          if ii==-1
            s = [ s ',' ];
          else
            s = [ s '] ' ];
          end
        end
        outstr.str(i,:) = s;
      end
    else                   % real intervals
      outstr.str = char(zeros(m,2*len+4));
      for i=1:m
        s = [ '[' ...
              dble2str_rnd(full(cinf(i)),commonexp,len,prec,expon,-1*InnerInclusion) ',' ...
              dble2str_rnd(full(csup(i)),commonexp,len,prec,expon,+1*InnerInclusion) '] ' ];
        outstr.str(i,:) = s;
      end
    end
    if rndold
      setround(rndold)
    end
    return
  end

  if nargout
    outstr = [];
    if c.complex
      for i=1:prod(size(cinf))
        outstr = [ outstr '[' ];
        for ii=[-1 1]
          if ii==-1
            x = full(cinf(i));
          else
            x = full(csup(i));
          end
          outstr = [ outstr dble2str_rnd(real(x),commonexp,len,prec,expon,ii) ];
          if imag(x)~=0
            strim = dble2str_rnd(imag(x),commonexp,len-1,prec,expon,ii);
            strim(isspace(strim)) = '';      % eliminate spaces
            if strim(1)~='-'
              strim = [ '+' strim ];
            end
            outstr = [ outstr strim 'i' ];
          end
          if ii==-1
            outstr = [ outstr ',' ];
          else
            outstr = [ outstr '] ' ];
          end
        end
      end
    else
      for i=1:prod(size(cinf))
        outstr = [ outstr '[' ...
                   dble2str_rnd(full(cinf(i)),commonexp,len,prec,expon,-1*InnerInclusion) ',' ...
                   dble2str_rnd(full(csup(i)),commonexp,len,prec,expon,+1*InnerInclusion) '] ' ];
      end
    end
    if rndold
      setround(rndold)
    end
    return
  end

  if issparse(c)
    [I,J] = find(spones(cinf)+spones(csup));    % inf or sup maybe zero
    if length(I)==0
      mid(c)
      if rndold
        setround(rndold)
      end
      return
    end
    if c.complex
      for i=1:length(I)
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        str = [ str '[' ];
        for ii=[-1 1]
          if ii==-1
            x = full(cinf(I(i),J(i)));
          else
            x = full(csup(I(i),J(i)));
          end
          strre = dble2str_rnd(real(x),commonexp,len,prec,expon,ii);
          if imag(x)==0
            chsign = ' ';
            strim = blanks(len);
          else
            signim = sign(imag(x));
            chsign = '+';
            if signim==-1
              chsign = '-';
            end
            strim = [ dble2str_rnd(abs(imag(x)), ...
                      commonexp,len-1,prec,expon,ii*signim) 'i' ];
          end
          str = [ str strre ' ' chsign strim ];
          if ii==-1
            str = [ str ',' ];
          else
            str = [ str '] ' ];
          end
        end
        disp( str )
      end
    else
      for i=1:length(I)
        str = sprintf('  (%d,%d)',I(i),J(i));
        str = [ str blanks(20-length(str)) ];
        str = [ str '[' ...
                   dble2str_rnd(full(cinf(I(i),J(i))),commonexp,len,prec,expon,-1*InnerInclusion) ',' ...
                   dble2str_rnd(full(csup(I(i),J(i))),commonexp,len,prec,expon,+1*InnerInclusion) '] ' ];
        disp( str )
      end
    end
    if rndold
      setround(rndold)
    end
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
          s = [ s '[' ];
          for ii=[-1 1]
            if ii==-1
              x = full(cinf(i,j));
            else
              x = full(csup(i,j));
            end
            strre = dble2str_rnd(real(x),commonexp,len,prec,expon,ii);
            if imag(x)==0
              chsign = ' ';
              strim = blanks(len);
            else
              signim = sign(imag(x));
              chsign = '+';
              if signim==-1
                chsign = '-';
              end
              strim = [ dble2str_rnd(abs(imag(x)), ...
                        commonexp,len-1,prec,expon,ii*signim) 'i' ];
            end
            s = [ s strre ' ' chsign strim ];
            if ii==-1
              s = [ s ',' ];
            else
              s = [ s '] ' ];
            end
          end
        end
        disp(s)
      end
    else                   % real intervals
      for i=1:m
        s = '';
        for j = j1:j2
          s = [ s '[' ...
                dble2str_rnd(full(cinf(i,j)),commonexp,len,prec,expon,-1*InnerInclusion) ',' ...
                dble2str_rnd(full(csup(i,j)),commonexp,len,prec,expon,+1*InnerInclusion) '] ' ];
        end
        disp(s)
      end
    end
    if loose & ~restricted, disp(' '); end
  end

  if rndold
    setround(rndold)
  end
