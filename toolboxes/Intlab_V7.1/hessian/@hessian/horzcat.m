function c = horzcat(varargin)
%HORZCAT      Implements  [a(1) a(2) ...]  for hessians
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = hessian(varargin{1});
  c.x = a.x;
  c.dx = a.dx;
  c.hx = a.hx;

  for i=2:length(varargin)
    a = hessian(varargin{i});
    c.x = [ c.x a.x ];
    c.dx = [ c.dx a.dx ];          % arrays stored columnwise
    c.hx = [ c.hx a.hx ];
  end

  c = class(c,'hessian');
