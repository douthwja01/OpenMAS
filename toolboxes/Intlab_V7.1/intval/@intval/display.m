function display(c,name,restricted)
%DISPLAY      Command window display of interval (rigorous)
%
%The three possibilities to display intervals may be defined to be default,
%  e.g. for x=midrad(3.14,.009); y=midrad(-5.333+.003i,1e-4); it is
%
%    » intvalinit('Display_'), x,y
%    intval x =
%        3.14__
%    intval y =
%      -5.333_ +  0.003_i
%
%    » intvalinit('DisplayInfsup'), x,y
%    intval x =
%    [    3.1309,    3.1491]
%    intval y =
%    [  -5.3332 +  0.0028i,  -5.3328 +  0.0032i]
%
%    » intvalinit('DisplayMidrad'), x,y
%    intval x =
%    <    3.1400,   0.0091>
%    intval y =
%    <  -5.3330 +  0.0030i,  0.0002>
%
%For details of the individual routines,
%
%  help disp_
%  help infsup
%  help midrad
%

%for internal use:
%  name                name of output variable
%  restricted == 1     no header, no extra lines output
%
%Call only with 1 or 3 input arguments
%

% written  11/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 08/26/12     S.M. Rump  global variables removed
%

  if nargin==1
    name = inputname(1);
    restricted = 0;
  elseif nargin==2
%     error('invalid call of display')
    restricted = 0;
  end

  INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
  if isequal(INTLAB_INTVAL_DISPLAY,'DisplayInfsup')
    infsup(c,name,restricted)
  elseif isequal(INTLAB_INTVAL_DISPLAY,'DisplayMidrad')
    midrad(c,name,restricted)
  else
    disp_(c,name,restricted)	    % covers disp_ and disp__
  end

  