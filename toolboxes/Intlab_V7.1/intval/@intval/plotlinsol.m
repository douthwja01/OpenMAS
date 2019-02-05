function Sigma = plotlinsol(A,b,printaxes,printvertexsol)
%PLOTLINSOL   Plots solution set of real interval linear system in 2 or 3 unknowns
%
%  Sigma = plotlinsol(A,b)
%
%CAREFUL: I had to use some dirty tricks to outsmart bugs in Matlab's
%convex hull routines. Please use only for thick input data (otherwise
%not much can be seen anyway).
%
%Plots the solution complex of a real interval linear system in 2 or 3 unknowns.
%Output Sigma is the interval hull of the solution complex. Output is NaN for
%singular input matrix.
%
%By specifying additional (optional) parameters, the x- and y-axes and/or all
%vertex solutions can be plotted, respectively:
%
%  plotlinsol(A,b,printaxes,printvertexsol)
%
%Here the vertex solutions is the set of all solutions of linear systems where
%the matrix and right hand side are vertices of the input interval system. If 
%all entries are thick, these are 2^(n^2)*2^n points. They are marked by
%asterisks. Careful, a little slow in 3 dimensions with 4096 points.
%
%Examples are an arrow-shaped solution set, here plotted with vertex solutions:
%
%  A = [infsup(2,4) infsup(-1,1);infsup(1,1) infsup(2,4)];  b = [infsup(-2,3);1];  plotlinsol(A,b,[],1)
%
%The solution set is convex in every orthant. In the following example this not immediately clear:
%
%  A = [infsup(2,4) infsup(-1,0);infsup(0,1) infsup(1,2)]; b = [infsup(-2,5);1]; plotlinsol(A,b)
%
%However, when plotting the axes convexity in every orthant becomes clear:
%
%  plotlinsol(A,b,1)
%
%Finally a modified batman:
%
%  A = [infsup(2,4) infsup(-1,1);infsup(-1,1) infsup(2,4)]; b = [infsup(-3,3);.8];  plotlinsol(A,b)
%
%Here is a 3-d example (courtesy of Walter Krämer, Wuppertal)
%
%  A = infsup([4 -2 -2;-1 5 -2;-1 -1 6],[5 2 1;2 5 2;2 2 6])
%  b = infsup([-2 -1 -.5]',[.5 1 2]')
%  plotlinsol(A,b)
%
%Finally, the cover page from Arnold's book:
%
%  A = ones(3)*infsup(0,2); A(1:4:end) = 3.5
%  b = ones(3,1)*infsup(-1,1)
%  plotlinsol(A,b)
%

% written  12/28/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 11/05/07     S.M. Rump  adapt to Matlab version 2007a
% modified 10/25/08     S.M. Rump  Matlab convhull bug (thanks to S. Shary for pointing to that)
% modified 10/10/09     S.M. Rump  3D and complete redesign
% modified 03/16/10     S.M. Rump  various Matlab bugs in convhull
%

  e = 1e-30;
  if 1+e==1-e                 % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  wng = warning;              % store warning mode    
  warning off

  A = intval(A);
  b = intval(b);

  if A.complex | b.complex
    error('input data must be real')
  end
    
  n = dim(A);
  if ( n~=2 ) & ( n~=3 )
    error('dimension must be 2 or 3')
  end
  
  if ~exist('printaxes') | isempty(printaxes)
    printaxes = 0;
  end
  
  if ~exist('printvertexsol') | isempty(printvertexsol)
    printvertexsol = 0;
  end
  
  if n==3
    Sigma = plotlinsol3(A,b,printaxes,printvertexsol);
  else
    midA = mid(A);
    radA = rad(A);
    midb = mid(b);
    radb = rad(b);

    % check non-singularity of interval matrix and convex hull of solution complex
    % following Theorem 5.1 (C4) in J. Rohn: Systems of Linear Interval Equations, LAA 126:39-78 (1989)
    wngold = warning;
    warning off
    P = [];
    for i=0:2^n-1
      S1 = diag((-1).^(dec2bin(i,n)-48));
      for j=0:2^n-1
        S2 = diag((-1).^(dec2bin(j,n)-48));
        d = diag(midA/(midA-S1*radA*S2));
        if any(d<=0.5)
          Sigma = NaN;
          if rndold, setround(rndold); end
          return
        end
        x = ( midA-S1*radA*S2 ) \ ( midb+S1*radb );
        P = [P x];
      end
    end
    warning(wngold)
    try
      index = convhull(P(1,:),P(2,:));
    catch
      %%% VVV correct Matlab bug: convhull doesn't like straight lines
      %%% dirty solution: add a little relative and absolute perturbation
      PP = P.*(1+1e-4*randn(size(P))) + 1e-4*randn(size(P))*max(abs(P(:)));
      index = convhull(PP(1,:),PP(2,:));
      %%% AAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    P = P(:,index);
    Xinf = [min(P(1,:));min(P(2,:))];
    Xsup = [max(P(1,:));max(P(2,:))];
    X = infsup(Xinf,Xsup);

    hold off
    fillbox(X,'r')
    v = axis;
    dx = v(2)-v(1);
    dy = v(4)-v(3);
    f = 0.01;
    axis(v + f*[-dx dx -dy dy])       % widen axes a little bit
    axis equal
    hold on

    orth = [ 0 0 ; 0 1 ; 1 0 ; 1 1 ];

    for k=1:2^n
      v = - (-1).^orth(k,:)';                 % vertex corresponding to orthant
      S = diag(v);                            % signature matrix corresponding to orthant
      [e,Y] = emptyintersect( X , hull(0,v*inf) );     % intersection of X with current orthant
      if ~any(e)
        exclude( [ midA-radA*S ; -midA-radA*S ] , [ -b.sup ; b.inf ] , Y );
      end
    end

    % make sure initial inclusion is dashed
    plotboxedges(X,'Color','w')
    plotboxedges(X,'LineStyle','--','Color','b')
    Sigma = X;

    if printaxes
      % AAAAAA  Matlab bug: X1 = X(1); X2 = X(2);
      T.type = '()';
      T.subs = {1};
      X1 = subsref(X,T);
      T.subs = {2};
      X2 = subsref(X,T);
      % VVVVVV Matlab bug: X1 = X(1); X2 = X(2);
      if in(0,X1)
        line([0 0],[X.inf(2);X.sup(2)],'Color','k','LineStyle',':');   % x-axis
      end
      if in(0,X2)
        line([X.inf(1);X.sup(1)],[0 0],'Color','k','LineStyle',':');   % y-axis
      end
    end

    if printvertexsol
      Ainf = A.inf;
      Asup = A.sup;
      binf = b.inf;
      bsup = b.sup;
      for i=0:2^(n*n)-1
        v = find(dec2bin(i,n*n) - 48);      % n*n-vector of 0/1, double('0') = 48
        a = Ainf;
        a(v) = Asup(v);
        for j=0:2^n-1
          w = find(dec2bin(j,n) - 48);      % n-vector of 0/1, double('0') = 48
          bb = binf;
          bb(w) = bsup(w);
          x = a\bb;
          if n==2
            plot(x(1),x(2),'k*')
          else
            plot3(x(1),x(2),x(3),'k*')
          end
        end
      end
    end
  end

  hold off
  warning(wng)          % restore warning mode
  if rndold
    setround(rndold)
  end
 
  
function Sigma = plotlinsol3(A,b,printaxes,printvertexsol) 
  n = dim(A);
  mA = mid(A);
  rA = rad(A);
  mb = mid(b);
  rb = rad(b);
  e = 1e-6*norm(mag(b),1);
  
  K = nchoosek(1:3*n,n);
  X = {};
  L = {};
  Mmin = zeros(1,n);
  Mmax = zeros(1,n);
  Sigma = intval(mA\mb);
  D = det(mA);
  for j=1:2^n
    XX = [];
    v = dec2bin(j,n) - 48;           % double('0') = 48
    ll = length(v);
    S = diag(2*v(end-n+1:end)-1);
    if ( D*det(mA-rA*S)<=0 ) | ( D*det(mA+rA*S)<=0 )
      Sigma = NaN;
      return
    end
    M = [ S ; rA*S-mA ; mA+rA*S ];
    B = [ zeros(n,1) ; -b.sup ; b.inf ];
    for i=1:size(K,1)
      x = M(K(i,:),:)\B(K(i,:));
%       if any(isnan(x)) | any(isinf(x)) 
%         Sigma = NaN;
%         if rndold, setround(rndold); end
%         return
%       end
      if ( M*x>=B-e ) & ( ~any(isnan(x)) ) & ( ~any(isinf(x)) ) 
        XX = [ XX ; x' ];
        Sigma = hull(Sigma,x);
      end
    end
    if ~isempty(XX)
      try
        if n==2
          L{end+1} = convhull(XX(:,1),XX(:,2));
        else
          L{end+1} = convhulln(XX);
        end
      catch
          %%% VVV correct Matlab bug: convhull doesn't like straight lines
          %%% dirty solution: add a little relative and absolute perturbation
          XX = XX.*(1+1e-4*randn(size(XX))) + 1e-4*randn(size(XX))*max(abs(XX(:)));
          L{end+1} = convhulln(XX);
          %%% AAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      X{end+1} = XX;
      Mmin = min(Mmin,min(XX));
      Mmax = max(Mmax,max(XX));
    end
  end
  
  close
  Mmin = Mmin*1.1 - .1*norm(Mmin);
  Mmax = Mmax*1.1 + .1*norm(Mmax);
  MM = [Mmin;Mmax];
  axis(MM(:)');
  axis equal
  hold on
  
  if printaxes
    h = line(1.3*[Mmin(1) 0 0;Mmax(1) 0 0],[0 Mmin(2) 0;0 Mmax(2) 0],[0 0 Mmin(3);0 0 Mmax(3)]);
    set(h,'Color','r')
  end
  
  for i=1:length(L)
    trisurf(L{i},X{i}(:,1),X{i}(:,2),X{i}(:,3));
  end
  
  if printvertexsol
    for i=1:2^(n^2)
      v = dec2bin(i,n^2) - 48;         % double('0') = 48
      ll = length(v);
      a = mA + reshape(2*v(end-n^2+1:end)-1,n,n).*rA;
      for j=1:2^n
        v = dec2bin(j,n) - 48;         % double('0') = 48
        ll = length(v);
        bb = mb + reshape(2*v(end-n+1:end)-1,n,1).*rb;
        x = a\bb;
        if n==2
          plot(x(1),x(2),'k*')
        else
          plot3(x(1),x(2),x(3),'k*')
        end
      end
    end
  end
  
  
function exclude(A,b,X)
% exclude all points in X with not(Ax+b<=0)

  if any( A*X+b > 0 )   % no feasible point
  %  fill(X.inf,X.sup,'w','EdgeColor','none')
    fillbox(X,'w','EdgeColor','none')
  else                                  % there are feasible points
    for i=1:size(A,1)
      exclude2(A(i,:),b(i),X.inf,X.sup)
    end
  end
  
  
function exclude2(a,b,Xinf,Xsup)
% exclude all points in [Xinf,Xsup] with not(ax+b<=0)

  P = [];
  XX = [Xinf(1) Xsup(1) Xsup(1) Xinf(1);Xinf(2) Xinf(2) Xsup(2) Xsup(2)];
  feasible = ( a*XX+b <= 0 );
  edges = [1 2;2 3;3 4;4 1];
  common = [2 1 2 1];
  uncommon = [1 2 1 2];
  for i=1:4
    I = edges(i,:);
    P = [P treatedge(a,b,XX(:,I),feasible(I),common(i),uncommon(i))];
  end
  if ~isempty(P)   %%% Matlab bug: fill may not work for points on a line
    try
      index = convhull(P(1,:),P(2,:));
      fill(P(1,index),P(2,index),'w','EdgeColor','none')
    end
  end
  
  
function P = treatedge(a,b,X,feasible,common,uncommon)
% treat edge X(:,1)..X(:,2)
  
  if feasible(1)==feasible(2)
    if feasible(1)==1                   % both feasible
      P = [];
    else                                % both not feasible
      P = X;
    end
  else                                  % one feasible, one not
    x = (-b-a(common)*X(common))/a(uncommon);
    Y = X(:,1);
    Y(uncommon) = x;
    P = [ Y X(:,find(feasible==0)) ];
  end
    
  
  
function fillbox(X,varargin)
% fill box [Xinf,Xsup] into figure subject to parameters in varargin

  Xinf = X.inf;
  Xsup = X.sup;
  fill([Xinf(1);Xsup(1);Xsup(1);Xinf(1)],[Xinf(2);Xinf(2);Xsup(2);Xsup(2)],varargin{:})

  
  
function plotboxedges(X,varargin)
% plot edges of box [Xinf,Xsup] into figure subject to parameters in varargin

  Xinf = X.inf;
  Xsup = X.sup;
  line([Xinf(1);Xsup(1)],[Xinf(2);Xinf(2)],varargin{:})
  line([Xsup(1);Xsup(1)],[Xinf(2);Xsup(2)],varargin{:})
  line([Xinf(1);Xsup(1)],[Xsup(2);Xsup(2)],varargin{:})
  line([Xinf(1);Xinf(1)],[Xinf(2);Xsup(2)],varargin{:})
