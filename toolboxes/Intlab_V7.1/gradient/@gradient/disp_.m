function disp_(c)
%DISP_        Display of interval gradients in "_" notation
%
%   disp_(c)
%

% written  06/04/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/05/12     S.M. Rump  complete redesign
%

  INTLAB_GRADIENT_NUMVAR = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  numvar = size(c.dx,2);
  if numvar~=INTLAB_GRADIENT_NUMVAR
    warning('**** number of dependent variables and partial derivatives do not coincide')
  end

  name = inputname(1);
  if isempty(name)                    % happens for display(gradientinit(random))
    name = 'ans';
  end

  INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
  setappdata(0,'INTLAB_INTVAL_DISPLAY','Display_');
  display(c,name)
  setappdata(0,'INTLAB_INTVAL_DISPLAY',INTLAB_INTVAL_DISPLAY);
    