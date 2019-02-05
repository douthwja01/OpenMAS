function display(p,name)
%DISPLAY      Command window display of polynom
%

% written  08/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 08/26/12     S.M. Rump  global variables removed
%

  loose = strcmp(get(0,'FormatSpacing'),'loose');
  if nargin==1
    name = inputname(1);
  end

  INTLAB_POLYNOM_DISPLAY = getappdata(0,'INTLAB_POLYNOM_DISPLAY');
  if ( size(p.e,2)<=1 ) & INTLAB_POLYNOM_DISPLAY   % univariate polynomial
    if loose, disp(' '); end
    if isa(p.c,'intval')
      disp([ 'intval polynom ' name '[' p.v '] = ' ])
      if loose, disp(' '); end
      display(p.c,'',1);
    else
      disp([ 'polynom ' name '[' p.v ']  = ' ])
      if loose, disp(' '); end
      disp(p.c)
    end
    if loose, disp(' '); end
  else                   % sparse display
    if size(p.e,2)==1    % univariate polynomial with sparse display
      p.e = ( max(p.e):-1:0 )';   % p.e is vector
      p.c = p.c(:);
      index = ( p.c==0 );
      p.e(index) = [];
      p.c(index) = [];
      if isempty(p.e)    % restore zero polynomial
        p.e = 0;
        p.c = typeadj( 0 , typeof(p.c) );
      end
      p.v = { p.v };
    end
    if loose, disp(' '); end
    if isa(p.c,'intval')
      INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
      s = [ 'intval polynom ' name '[' p.v{1} ];
      for i=2:size(p.e,2)
        s = [ s ',' p.v{i} ];
      end
    else
      s = [ 'polynom ' name '[' p.v{1} ];
      for i=2:size(p.e,2)
        s = [ s ',' p.v{i} ];
      end
    end
    disp([ s '] = ' ])
    [ p.e index ] = sortrows(p.e);
    p.e = flipud(p.e);
    p.c = flipud(p.c(index));
    if loose, disp(' '); end
    out = disp2str(p.c);
    disp(out.exp)
    pe_max = max(p.e,[],1);
    pe_max(pe_max==0) = 1;
    k = floor(log10(pe_max)) + 1;
    for j=1:size(p.e,2)
      format{j} = [ '%' int2str(k(j)) 'd' ];
    end
    for i=1:size(p.e,1)
      s = [ out.str(i,:) '  ' ];
      for j=1:size(p.e,2)
        if p.e(i,j)>1
          s = [ s p.v{j} '^' sprintf(format{j},p.e(i,j)) '  ' ];
        elseif p.e(i,j)==1
          s = [ s p.v{j} blanks(k(j)+3) ];  
        elseif p.e(i,j)==0
          s = [ s blanks(length(p.v{j})+k(j)+3) ];
        end
      end
      disp(s)
    end
    if loose, disp(' '); end
  end

  