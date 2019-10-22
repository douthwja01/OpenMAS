function c = vertcat(varargin)
%VERTCAT      Implements  [a(1) ; a(2) ; ...]  for hessians
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = hessian(varargin{1});
  c.x = a.x.';
  [m n] = size(a.x);
  index = reshape( 1:(m*n) , m , n )';
  c.dx = a.dx( :, index(:) );
  c.hx = a.hx( :, index(:) );

  for i=2:length(varargin)
    a = hessian(varargin{i});
    c.x = [ c.x a.x.' ];
    [m n] = size(a.x);
    index = reshape( 1:(m*n) , m , n )';
    c.dx = [ c.dx a.dx( : , index(:) ) ];       % arrays stored columnwise
    c.hx = [ c.hx a.hx( : , index(:) ) ];
  end
  [m n] = size(c.x);
  index = reshape( 1:(m*n) , m , n )';
  c.x = c.x.';
  c.dx = c.dx( : , index(:) );
  c.hx = c.hx( : , index(:) );

  c = class(c,'hessian');
