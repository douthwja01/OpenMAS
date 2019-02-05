function display(c,name)
%DISPLAY      Command window display of gradient
%
%Gradients of row vectors are 3-dimensional arrays. Therefore the gradient .dx
%  of a sparse array with more than one column is displayed as full gradient.
%

%Second parameter name for internal purposes
%

% written  10/16/98     S.M. Rump
% modified 11/06/99     S.M. Rump  omit ans
% modified 03/07/04     S.M. Rump  output changed to gradients of individual components
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/05/12     S.M. Rump  comment to sparse data
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_GRADIENT_NUMVAR = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  if nargin<2
    name = inputname(1);
    if isempty(name)                    % happens for display(gradientinit(random))
      name = 'ans';
    end
  end

  numvar = size(c.dx,2);
  if numvar~=INTLAB_GRADIENT_NUMVAR
    warning('**** number of dependent variables and partial derivatives do not coincide')
  end

  sizecx = size(c.x);
  if isa(c.x,'intval')
    strintval = 'intval ';
  else
    strintval = '';
  end
    
  % display .x
  if loose, disp(' '); end
  disp([ strintval 'gradient value ' name '.x = ' ])
  if loose, disp(' '); end
  display_(c.x,strintval,name)
  
  % display .dx
  if loose, disp(' '); end
  disp([ strintval 'gradient derivative(s) ' name '.dx = ' ])
  if ( sizecx(2)==1 ) & ( length(sizecx)==2 )   % scalar or column vector
    display_(c.dx,strintval,name);
  else
    index = ones(1,length(sizecx));
    for i=1:prod(sizecx)
      temp = reshape(c.dx(i,:),1,numvar);
      str = [strintval name '.dx('];
      for j=1:length(index)
        str = [str int2str(index(j)) ','];
      end
      disp([str ':) = '])
      display_(temp,strintval,name);
      index = nextindex(index,sizecx);
    end
  end
  
  if loose, disp(' '); end
  
  if rndold
    setround(rndold)
  end

  
function display_(c,strintval,name)
% display dependent on strintval
  if isempty(strintval)
    % eval([name ' = c;']);       % make sure name is displayed
    % eval(['disp(' name ')'])
    disp(c)
  else
    display( c , name , 1 )
  end
    
  
function index = nextindex(index,size)
% compute next multi-index subject to size
  i=1;
  while i<=length(size)
    index(i) = index(i)+1;
    if index(i)<=size(i)
      return
    end
    index(i) = 1;
    i = i+1;
  end
  