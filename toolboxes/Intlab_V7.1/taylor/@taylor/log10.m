function a = log10(a)
%LOG10        Taylor logarithm  log10(a)
%

% written  05/22/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  a = log(a);
  if isa(a.t,'intval')
    INTLAB_STDFCTS_LOG10_ = getappdata(0,'INTLAB_STDFCTS_LOG10_');
    a.t = a.t .* infsup(INTLAB_STDFCTS_LOG10_.INF,INTLAB_STDFCTS_LOG10_.SUP);
  else
    a.t = a.t / log(10);
  end
  