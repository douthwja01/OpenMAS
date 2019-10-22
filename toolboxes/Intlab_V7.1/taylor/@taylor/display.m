function display(c,name)
%DISPLAY      Command window display of Taylor
%

%Second parameter name for internal purposes
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/07/12     S.M. Rump  complete redesign
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_TAYLOR_ORDER = getappdata(0,'INTLAB_TAYLOR_ORDER');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  if nargin<2
    name = inputname(1);
    if isempty(name)                    % happens for display(taylorinit(random))
      name = 'ans';
    end
  end

  numvar = size(c.t,1)-1;
  if numvar~=INTLAB_TAYLOR_ORDER
    warning('**** number of dependent variables and partial derivatives do not coincide')
  end

  if isa(c.t,'intval')
    strintval = 'intval ';
  else
    strintval = '';
  end
      
  if ( c.size(2)==1 ) & ( length(c.size)==2 )   % scalar or column vector
    display_(c.t.',strintval,name);
  else
    NN = size(c.t,1);
    index = ones(1,length(c.size));
    for i=1:prod(c.size)
      temp = reshape(c.t(:,i),1,NN);
      str = [strintval name '.t('];
      for j=1:length(index)
        str = [str int2str(index(j)) ','];
      end
      disp([str ':) = '])
      display_(temp,strintval,name);
      index = nextindex(index,c.size);
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
  