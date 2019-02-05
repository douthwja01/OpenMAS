function out = disp2str(x)
%DISP2STR     Output of intval x into out
%
%The output produced by disp_(x(:)), infsup(x(:)), or midrad(x(:)),
%  according to the format in use, is returned into  out  as follows:
%
% out.exp  empty if no common exponent printed, otherwise
%            string of common exponent
% out.str  column array of strings representing to x(:)
%

% written  08/29/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_INTVAL_DISPLAY = getappdata(0,'INTLAB_INTVAL_DISPLAY');
%VVVV x = x(:);
  x = reshape(x,prod(size(x)),1);
%AAAA Matlab V5.2 bug fix
  switch INTLAB_INTVAL_DISPLAY
    case 'DisplayInfsup', out = infsup(x,[],[]);
    case 'DisplayMidrad', out = midrad(x,[],[]);
    case 'Display_', out = disp_(x,[],[]);
  end
