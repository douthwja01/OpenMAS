function plot(varargin)
%PLOT         Plot X vs. Y for real interval vectors X or Y
%
%Call
%
%  plot(Y)   or   plot(X,Y)   or   plot(X,Y,c)
%
%For Y being a real interval vector, the curves (x,Y.inf) and (x,Y.sup) are filled and plotted.
%For X being a real interval vector, the boxes X(i) x Y(i) are plotted.
%
%The area between curves or boxes is filled with color c; default is 'r' for red.
%

% written  11/15/07     S.M. Rump  (inspired by R. Malti)
%

  if nargin==1
    c = 'r';
    Y = varargin{1};
    if size(Y,1)==1
      Y = Y.';
    end
    if isreal(Y)
      X = (1:size(Y,1))';
    else                     % plot real against imaginary part
      X = real(Y);
      Y = imag(Y);
    end
  elseif nargin==2
    X = varargin{1};
    Y = varargin{2};
    if ischar(Y)             % input X not specified
      c = Y;
      if size(X,1)==1
        X = X.';
      end
      if isreal(X)
        Y = X;
        X = (1:size(X,1))';
      else                   % plot real against imaginary part
        Y = imag(X);
        X = real(X);
      end
    else
      c = 'r';
    end
  else                      % input X,Y,c
    X = varargin{1};
    Y = varargin{2};
    c = varargin{3};    
  end
  if size(X,1)==1
    X = X';
    Y = Y';
  end

  % prepare axis
  if isa(X,'intval')
    xmin = min(X.inf(:));
    xmax = max(X.sup(:));
  else
    xmin = min(X(:));
    xmax = max(X(:));
  end
  if isa(Y,'intval')
    ymin = min(Y.inf(:));
    ymax = max(Y.sup(:));
  else
    ymin = min(Y(:));
    ymax = max(Y(:));
  end
  % take care of zero width
  dx = (xmax-xmin)*0.05;
  xmin = xmin - dx - realmin;
  xmax = xmax + dx + realmin;
  dy = (ymax-ymin)*0.05;
  ymin = ymin - dy - realmin;
  ymax = ymax + dy + realmin;
  axis([xmin xmax ymin ymax])
  plot(xmin+.5*(xmax-xmin),ymin+.5*(ymax-ymin),'w');
  hold on

  if isa(X,'intval')
    if X.complex
      error('Plot for intervals only for real input data.')
    end
    if isa(Y,'intval')      % interval vectors X and Y
      if Y.complex
        error('Plot for intervals only for real input data.')
      end
      X.inf = X.inf(:)';
      X.sup = X.sup(:)';
      Y.inf = Y.inf(:)';
      Y.sup = Y.sup(:)';
      fill([X.inf;X.sup;X.sup;X.inf],[Y.inf;Y.inf;Y.sup;Y.sup],c)
    else                    % interval vector X, point vector Y
      X.inf = X.inf(:)';
      X.sup = X.sup(:)';
      Y = Y(:)';
      fill([X.inf;X.sup],[Y;Y],c)
%       fill([X.inf;flipud(X.sup)],[Y;flipud(Y)],c)
    end
  else                      % point vector X, interval vector Y
    if Y.complex
      error('Plot for intervals only for real input data.')
    end
    fill([X;flipud(X)],[Y.inf;flipud(Y.sup)],c)
  end
  hold off
  