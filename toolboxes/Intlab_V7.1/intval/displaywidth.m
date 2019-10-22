function res = displaywidth(digits,see)
%DISPLAYWIDTH Initialization of width of display
%
%   digits = displaywidth     gets current width of display, default 120
%   displaywidth(digits)      sets width of display to width, minimum 110
%
%When changing the current value, a corresponding message will be printed.
%To suppress the message, use
%
%   digits = displaywidth(digits,0)
%

% written  12/15/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_DISPLAY_WIDTH = getappdata(0,'INTLAB_DISPLAY_WIDTH');
  if ( nargin==0 )                  % getting current width
    res = INTLAB_DISPLAY_WIDTH;
  else                              % setting current width
    INTLAB_DISPLAY_WIDTH = max(digits,110);
    setappdata(0,'INTLAB_DISPLAY_WIDTH',INTLAB_DISPLAY_WIDTH);
    res = INTLAB_DISPLAY_WIDTH;
    if nargin==1
      see = 1;
    end
    if see
      disp(['===> width of display ' int2str(res) ' characters'])
    end
  end
