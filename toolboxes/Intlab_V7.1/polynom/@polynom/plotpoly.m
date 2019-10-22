function [X,y,z] = plotpoly(p,xl,xr,kmax)
%PLOTPOLY     Plot routine for real (interval) polynomial in one or two variables
%
%For univariate polynomial p, 
%   plotpoly(p)
%plots p within [-r,r], where r is a rootbound for p. For interval polynomial p,
%the lower and upper bound is plotted in blue and red, respectively.
%   plotpoly(p,xl,xr)
%plots p within [xl,xr], both with 100 meshpoints. 
%   [x,y] = plotpoly(p,xl,xr,kmax)
%plots p within [xl,xr] with kmax meshpoints. Optional output parameters store the
%vectors of x- and y-values, where y is an interval vector for interval polynomial p.
%The x-axis is displayed if within the plot.
%
%For polynomial p in two unknowns, 
%   plotpoly(p,lb,ub,kmax)
%plots p within lb and ub, where lb and ub are two-vectors specifying the lower and upper
%bound for the two variables, respectively. Input parameter kmax is optional, default is 100.
%Input kmax may also be two-vector specifiying number of meshpoints for first and second variable.
%   Z = plotpoly(p,lb,ub)
%gives back matrix Z of polynomial values, and 
%   [X,Y,Z] = plotpoly(p,lb,ub)
%gives back vectors X and Y of values of first and second variable, respectively, and matrix Z of
%polynomial values, such that surf(X,Y,Z) produces polynomial plot (for interval polynomial use
%  surf(X,Y,Z.inf), hold on, surf(X,Y,Z.sup), hold off ).
%

% written  09/14/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 11/05/07     S.M. Rump  zero line added
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isreal(p)
    error('polynomial must be real (point or interval)')
  end
  
  if size(p.e,2)==1                   % univariate polynomial
    
    if nargin==1
      xl = -rootbound(p);
      xr = -xl;
    end
    
    if isnan(xl) | ( xr-xl==0 )
      xl = -1;
      xr = 1;
    end
    
    if nargin<4
      kmax = 100;
    end
    
    x = linspace(xl,xr,kmax);
    if nargout>0
      X = x;
    end
    if isa(p.c,'intval')                 % univariate interval polynomial
      y = polyval(p,x);
      plot( x , y.inf , 'b', x , y.sup , 'r' )
    else                                 % non-interval polynomial
      plot( x , polyval(p.c,x) , 'r' )
    end
    
    A = axis;
    if ( A(3)<=0 ) & ( 0 <= A(4) )
      hold on
      plot(A(1:2),[0 0])
      hold off
    end
    
  elseif size(p.e,2)==2                  % multivariate polynomial in two unknowns
    
    if nargin<4
      xmax = 100;
      ymax = 100;
    else
      if length(kmax)==1
        xmax = kmax;
        ymax = kmax;
      elseif length(kmax)==2
        xmax = kmax(1);
        ymax = kmax(2);
      else
        error('invalid parameter kmax')
      end
    end
    if length(xl)==1
      xl = repmat(xl,1,2);
    end
    if length(xr)==1
      xr = repmat(xr,1,2);
    end
    if ( length(xl)~=2 ) | ( length(xr)~=2 )
      error('invalid input parameters for plotpoly')
    end
    x = linspace(xl(1),xr(1),xmax);
    y = linspace(xl(2),xr(2),ymax);
    xx = repmat(x,ymax,1);
    yy = repmat(y',1,xmax);
    z = reshape( polyval(p,[xx(:) yy(:)]) , ymax , xmax );
    
    if isa(p.c,'intval')
      surf(x,y,z.inf)
      hold on
      surf(x,y,z.sup)
      hold off
    else
      surf(x,y,z)
    end
    xlabel(p.v(1));
    ylabel(p.v(2));
    
    if nargout==1
      X = z;
    elseif nargout~=0
      X = x;
    end
    
  else
    error('polynomial plot only for polynomials in one or two unknowns')
  end
  
  setround(rndold)
