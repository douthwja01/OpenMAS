function y = log2(x)
%LOG2         Implements  log2(x)  for intervals (binary logarithm)
%
%   y = log2(x)
%
%interval standard function implementation
%

% written  11/08/07     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_LOG2_ = getappdata(0,'INTLAB_STDFCTS_LOG2_');

  % log(x) * ( 1/log(2) )
  y = log(x) * infsup(INTLAB_STDFCTS_LOG2_.INF , ...
                      INTLAB_STDFCTS_LOG2_.SUP );
  index = ( x==0 );           
  if any(index(:))
    %VVVV  y(index) = -inf;
    s.type = '()'; s.subs = {index}; y = subsasgn(y,s,-inf);
    %AAAA  Matlab bug fix    
  end
