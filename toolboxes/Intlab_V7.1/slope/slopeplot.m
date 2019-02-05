function slopeplot(f,xs,X,Xplot,n)
%SLOPEPLOT    Plot of one-dimensional function with slope expansion
%
%  slopeplot(f,xs,X,Xplot,n)
%
%Simple functions, which can be written in one string, can be entered
%  directly. The unknown must be 'x'. E.g.,
%
%  slopeplot('sqrt(abs(x))',2,infsup(-1,1));            or
%  slopeplot('sin(x)',0,infsup(-1,1));                  or
%  slopeplot('x.*x',1,infsup(-2,3));
%
%For a one-dimensional function f to be called by u=f(x) consisting of
%  operators overloaded by slope operators such that for a column vector
%  v, f(v) computes the vector f(v(i)), i=1:length(v).
%  Then the function f is expanded within X with respect to xs and
%  plotted together with the slopes in the interval Xplot).
%
%Parameter Xplot is optional; default is  hull(xs,X) * midrad(1,0.2);
%Parameter n is optional (number of grid points for function plot); default is 100
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if prod(size(xs))~=1
    error('slopeplot only for one-dimensional functions')
  end

  if nargin==3 | isempty(Xplot)
    Xplot = hull(xs,X);
    Xplot = Xplot + midrad(0,0.2*diam(Xplot));
  end

  if nargin<=4
    n = 100;
  end

  % slope evaluation
  x = slopeinit(xs,X);
  eval([ 'u = ' f ';' ])

  % plot function
  v = linspace(Xplot.inf,Xplot.sup,n);
  x = v;
  eval([ 'y = ' f ';' ])
  plot(x,y,'k')
  hold on

  % plot expansion range X
  A = axis;
  dy = A(4)-A(3);
  axis([ A(1:2) A(3)-.1*dy A(4)+.1*dy ])
  H = (A(4)-A(3))/200;         % height of box for X
  infX = inf(X);
  supX = sup(X);
  XX = [ infX supX supX infX ];
  patch(XX,A(3)+[-H -H H H],'b')

  % plot expansion point xs
  if isa(xs,'intval')
    infxs = inf(xs);
    supxs = sup(xs);
  else
    infxs = xs;
    supxs = xs;
  end
  xxs = [ infxs supxs supxs infxs ];
  patch(xxs,A(3)+2*[-H -H H H],'r')
  plot([infX infxs supxs],A(3)*ones(1,3),'b')

  % lines infX,supX,xs to function
  x = infX; eval(['y = ' f ';']); plot([infX infX],[A(3) y],'b:')
  x = supX; eval(['y = ' f ';']); plot([supX supX],[A(3) y],'b:')
  x = infxs; eval(['y = ' f ';']); plot([infxs infxs],[A(3) y],'r:')
  x = supxs; eval(['y = ' f ';']); plot([supxs supxs],[A(3) y],'r:')
  x = infxs; eval(['y = ' f ';']); plot(infxs,y,'r.')
  x = supxs; eval(['y = ' f ';']); plot(supxs,y,'r.')

  % slopes
  plot( v , inf(u.c) + inf(u.s)*(v-infxs) , 'b-' );
  plot( v , inf(u.c) + sup(u.s)*(v-infxs) , 'b-' );
  plot( v , sup(u.c) + inf(u.s)*(v-supxs) , 'b-' );
  plot( v , sup(u.c) + sup(u.s)*(v-supxs) , 'b-' );

  % add text X and xs
  if infxs+.5*(supxs-infxs)<infX+.5*(supX-infX)
    set(gcf,'DefaultTextColor','blue')
    text(supX,A(3)-8*H,'\bf X')
    set(gcf,'DefaultTextColor','red')
    text(infxs,A(3)-8*H,'\bf xs')
  else
    set(gcf,'DefaultTextColor','blue')
    text(infX,A(3)-8*H,'\bf X')
    set(gcf,'DefaultTextColor','red')
    text(supxs,A(3)-8*H,'\bf xs')
  end

  % add title
  set(gcf,'DefaultTextColor','default')
  strX = infsup(X);
  strX(isspace(strX)) = '';
  if isa(xs,'intval')
    strxs = infsup(xs);
  else
    strxs = sprintf('%g',xs);
  end
  strxs(isspace(strxs)) = '';
  title(['slope expansion of ' '' f '' ' in X=' strX ...
         ' with respect to xs=' strxs ])

  hold off
  
  setround(rndold)
