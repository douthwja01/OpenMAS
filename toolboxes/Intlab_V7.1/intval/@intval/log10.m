function y = log10(x)
%LOG10        Implements  log10(x)  for intervals (logarithm to base 10)
%
%   y = log10(x)
%
%interval standard function implementation
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/07/04     S.M. Rump  accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/03/05     S.M. Rump  completely simplified
% modified 10/18/08     S.M. Rump  out-of-range flag
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_LOG10_ = getappdata(0,'INTLAB_STDFCTS_LOG10_');

  % log(x) * ( 1/log(10) )
  y = log(x) * infsup(INTLAB_STDFCTS_LOG10_.INF , ...
                      INTLAB_STDFCTS_LOG10_.SUP );
  index = ( x==0 );           
  if any(index(:))
    %VVVV  y(index) = -inf;
    s.type = '()'; s.subs = {index}; y = subsasgn(y,s,-inf);
    %AAAA  Matlab bug fix    
  end
