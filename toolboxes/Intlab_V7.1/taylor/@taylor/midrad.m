function midrad(c)
%MIDRAD       Display of interval Taylor in midrad notation
%
%   midrad(c)
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/07/12     S.M. Rump  complete redesign
%

  INTLAB_TAYLOR_ORDER = getappdata(0,'INTLAB_TAYLOR_ORDER');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  name = inputname(1);
  if isempty(name)                    % happens for display(taylorinit(random))
    name = 'ans';
  end

  numvar = size(c.t,1)-1;
  if numvar~=INTLAB_TAYLOR_ORDER
    warning('**** number of dependent variables and partial derivatives do not coincide')
  end
  
  INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
  setappdata(0,'INTLAB_INTVAL_DISPLAY','DisplayMidrad');
  display(c,name)
  setappdata(0,'INTLAB_INTVAL_DISPLAY',INTLAB_INTVAL_DISPLAY);
