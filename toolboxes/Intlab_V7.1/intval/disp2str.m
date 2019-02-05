function out = disp2str(x)
%DISP2STR     Output of x into out
%
%The output produced by disp(x(:)), according to the format in use,
%  is returned into  out  as follows:
%
% out.exp  empty if no common exponent printed, otherwise
%            string of common exponent
% out.str  column array of strings representing to x(:)
%

% written  08/29/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  format = get(0,'format');
  if ~isequal(format,'short') & ~isequal(format,'long') & ...
     ~isequal(format,'shortE') & ~isequal(format,'longE')
    format = 'short';
  end

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');
  commonexp = 0;
  if ~isequal(format(end),'E')
    % calculate exponent range and common factor
    if isreal(x)
      xmax = abs(x);
      xmax(find(isinf(xmax))) = 0;              % fast for sparse matrices
      xmax = max(xmax(:));
    else
      xre = abs(real(x));
      xim = abs(imag(x));
      xre(find(isinf(xre))) = 0;                % fast for sparse matrices
      xim(find(isinf(xim))) = 0;
      xremax = max(xre(:));
      ximmax = max(xim(:));
      xmax = max(xremax,ximmax);
    end
    % care for sparse matrices
    xmax = full(xmax);
    if isempty(xmax)
      xmax = 0;
    end

    emax = floor(log10( xmax + (xmax==0) ));
    if isreal(x)
      if isequal(format,'short') & emax >= 3
        commonexp = emax;
      elseif isequal(format,'long') & emax >= 2
        commonexp = emax;
      end
      if emax <= -4
        commonexp = emax+1;
      end
    else
      if emax >= 2
        commonexp = emax;
      end
      if emax <= -3
        commonexp = emax+1;
      end
    end
  end
  x = x / INTLAB_INTVAL_POWER10.sup(1,commonexp+341);   % 10^commonexp

  if isreal(x)
    switch format
      case 'short',  len = 10; prec = 4;  expon = 0;
      case 'long',   len = 19; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 24; prec = 15; expon = 1;
    end
  else
    switch format
      case 'short',  len = 9;  prec = 4;  expon = 0;
      case 'long',   len = 18; prec = 14; expon = 0;
      case 'shortE', len = 13; prec = 4;  expon = 1;
      case 'longE',  len = 23; prec = 15; expon = 1;
    end
  end

  nblanks = 0;
  if commonexp~=0                      % common factor
    strexp = sprintf('%04d',commonexp);
    if commonexp >= 0
      strexp(1) = '+';
    end
    out.exp = [ '  1.0e' strexp ' *' ];
  else
    out.exp = [];
  end

  n = length(x(:));
  format = [ '%' sprintf('%d',len) '.' sprintf('%d',prec) 'e' ];
  if ~expon
    format(end) = 'f';
  end

  if isreal(x)            % real input
    for i=1:n
      out.str(i,:) = sprintf(format,x(i));
    end
  else                   % complex input
    for i=1:n
      chsign = ' +';
      if imag(x(i))<0
        chsign = ' -';
      end
      out.str(i,:) = [ sprintf(format,real(x(i))) chsign       ...
                       sprintf(format,abs(imag(x(i)))) 'i' ];
    end
  end
  
  if rndold
    setround(rndold)
  end
