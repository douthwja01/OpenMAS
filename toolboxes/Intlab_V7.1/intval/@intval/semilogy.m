function semilogy(varargin)
%SEMILOGY     Semilogarithmic plot X vs. Y for real interval vectors X or Y
%
%Call
%
%  semilogy(Y)   or   semilogy(X,Y)   or   semilogy(X,Y,c)
%
%The y-axis is in logarithmic scale.
%
%For Y being a real interval vector, the curves (x,Y.inf) and (x,Y.sup) are plotted.
%For X being a real interval vector, the boxes X(i) x Y(i) are plotted.
%
%The area between curves or boxes is filled with color c; default is 'r' for red.
%

% written  11/15/07     S.M. Rump  (inspired by R. Malti)
%

  plot(varargin{:})
  set(gca,'yscale', 'log')
