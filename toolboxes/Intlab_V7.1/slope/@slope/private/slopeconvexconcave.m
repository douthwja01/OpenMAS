function us = slopeconvexconcave(f,strfs,a,convex)
%SLOPECONVEXCONCAVE  Improved slope for convex and concave functions
%
%Internal function
%

%Input  f       string such that feval(f,a) evaluates f(a)
%       strfs   string such that eval(str) gives back f'(a) for
%                 str = strrep(strfs,'%','a')
%       a       input argument
%       convex  1/0 for convex/concave function
%
%Output us      improved slope
%

% written  12/06/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  indexc = 1:INTLAB_SLOPE.NUMVAR;
  indexr = 2:INTLAB_SLOPE.NUMVAR+1;
  xs = a.r(:,indexc);
  X = a.r(:,indexr);

  ws = warning;
  warning off

  if convex
    bndX = intval(inf(X));
    bndxs = intval(inf(xs));
  else
    bndX = intval(sup(X));
    bndxs = intval(sup(xs));
  end
  fbndX = feval(f,bndX);
  fbndxs = feval(f,bndxs);
  N = bndX - bndxs;
  s = ( fbndX - fbndxs ) ./ N ;
  index = ( N.inf<=0 ) & ( N.sup>=0 );
  if any(index(:))
    s(index) = infsup(-inf,inf);
  end
  hullXxs = hull(bndX,bndxs);
  strfs = strrep(strfs,'%','hullXxs');
  eval(['s1 = max( inf(s) , inf(' strfs ') );'])

  if convex
    bndX = intval(sup(X));
    bndxs = intval(sup(xs));
  else
    bndX = intval(inf(X));
    bndxs = intval(inf(xs));
  end
  fbndX = feval(f,bndX);
  fbndxs = feval(f,bndxs);
  N = bndX - bndxs;
  s = ( fbndX - fbndxs ) ./ N ;
  index = ( N.inf<=0 ) & ( N.sup>=0 );
  if any(index(:))
    s(index) = infsup(-inf,inf);
  end
  hullXxs = hull(bndX,bndxs);
% strfs = strrep(strfs,'%','hullXxs');   %% already defined %%
  eval(['s2 = min( sup(s) , sup(' strfs ') );'])

  warning(ws)
  us = infsup(s1,s2) .* a.s;
