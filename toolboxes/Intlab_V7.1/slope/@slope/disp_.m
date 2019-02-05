function disp_(u,str)
%DISP_        Display of slope variables in "_" notation
%
%  disp_(u)
%

%Second parameter name for internal purposes: display full structure
%

% written  06/04/09     S.M. Rump
% modified 02/28/10     S.M. Rump  vector output and multi-dimensional arrays
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  numvar = size(u.s,2);
  if numvar~=INTLAB_SLOPE.NUMVAR
    warning('**** number of dependent variables in slope expansion inconsitent')
  end

  if nargin==2
    if isequal(str,'internal')
      disp([ 'slope range (internal structure) ' inputname(1) '.r = ' ])
      display(u.r,'',1);
      disp([ 'slope slope ' inputname(1) '.s = ' ])
      display(u.s,'',1);
      return
    else
      error('Invalid call of slope disp_')
    end
  end

  if loose, disp(' '); end
  disp([ 'slope intval center ' inputname(1) '.c = ' ])
  if loose, disp(' '); end
  ur = u.r(:,1);
  disp_( reshape(ur,u.size) , '' , 1 )
  if loose, disp(' '); end

  if loose, disp(' '); end
  disp([ 'slope intval range ' inputname(1) '.r = ' ])
  if loose, disp(' '); end
  ur = u.r(:,INTLAB_SLOPE.NUMVAR+1);
  disp_( reshape(ur,u.size) , '' , 1 )
  if loose, disp(' '); end

  if loose, disp(' '); end
  disp([ 'slope intval slope ' inputname(1) '.s = ' ])
  if loose, disp(' '); end
  if ( length(u.size)==2 ) & ( u.size(2)==1 )
    disp_( reshape(u.s,[u.size(1) numvar]) , '' , 1 )
  else
    disp_( reshape(u.s,[u.size numvar]) , '' , 1 )
  end
  if loose, disp(' '); end

  