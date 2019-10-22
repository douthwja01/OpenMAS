function display_gen(c,name,display_routine)
%DISPLAY_GEN  General display routine for hessian (covers disp_, infsup, midrad)
%

%Second parameter name for internal purposes
%

% written  04/04/04     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_HESSIAN_NUMVAR = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  numvar = size(c.dx,1);   
  if numvar~=INTLAB_HESSIAN_NUMVAR      % this should never happen
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
  disp([ strintval 'Hessian value ' name '.x = ' ])
  if loose, disp(' '); end
  feval(display_routine,c.x,strintval,name)
  
  % display .dx
  if loose, disp(' '); end
  disp([ strintval 'Hessian first derivative(s) ' name '.dx = ' ])
  if ( sizecx(2)==1 ) & ( length(sizecx)==2 )   % scalar or column vector
    feval(display_routine,reshape(transpose(c.dx),[sizecx(1) numvar]),strintval,name);
  else
    index = ones(1,length(sizecx));
    for i=1:prod(sizecx)
      temp = reshape(c.dx(:,i),1,numvar);
      str = [strintval name '.dx('];
      for j=1:length(index)
        str = [str int2str(index(j)) ','];
      end
      disp([str ':) = '])
      feval(display_routine,temp,strintval,name);
      index = nextindex(index,sizecx);
    end
  end
  
  % display .hx
  if loose, disp(' '); end
  disp([ strintval 'Hessian second derivative(s) ' name '.hx = ' ])
  if loose, disp(' '); end
  nn = prod(sizecx);
  if nn==1                  % one variable
    temp = reshape(c.hx,numvar,numvar);
    temp = temp + transpose(temp);
    feval(display_routine,temp,strintval,name);
  else                      % several variables
    index = ones(1,length(sizecx));
    for i=1:nn
      % cures Matlab bug:  a=sparse([],[],[],1,1); reshape(a,1,1)*2  is not zero
      if any(find(c.hx(:,i)))
        temp = reshape(c.hx(:,i),numvar,numvar);
        temp = temp + transpose(temp);
      else
        if issparse(c.hx)
          temp = sparse([],[],[],numvar,numvar);
        else
          temp = zeros(numvar);
        end
        if isa(c.hx,'intval')
          temp = intval(temp);
        end
      end
      str = [strintval name '.hx('];
      for j=1:length(index)
        str = [str int2str(index(j)) ','];
      end
      disp([str ':,:) = '])
      feval(display_routine,temp,strintval,name);
      index = nextindex(index,sizecx);
    end
  end
  
  if loose, disp(' '); end

  
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
  