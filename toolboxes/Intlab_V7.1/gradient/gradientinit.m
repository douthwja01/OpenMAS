function r = gradientinit(v)
%GRADIENTINIT Initialization of gradient variable
%
%  x = gradientinit(v)
%
%The dependent variable(s) x are identified and initialized with value(s) v
%  v may be scalar for one independent variable, or vector (or matrix)
%  for several independent variables
%
%The call
%
%  gradientinit
%
%without input parameters sets the number of dependent variables to zero.
%Hessians are frequently sparsely populated. Therefore, one can choose to 
%store the derivative dense or sparse, see "help sparsegradient". 
%
%The gradient of arrays is stored in the 'next dimension'. So y.dx is 
%  3-dimensional for a gradient row vector y. Since Matlab does not 
%  support multi-dimensional sparse arrays, y.dx is not accessible for
%  sparse arrays with more than one column.
%
%As a simple example of gradients
%
%  u = gradientinit([ -3 ; 3+4i ])
%
%initializes the gradient package to have two dependent variables. The
%variable u, a column vector, with values u.x(1) = -3 and u.x(2) = 3+4i 
%has gradients u(1).dx = [1 0] and u(2).dx = [0 1].
%If after that the statement
%
%  v = gradient( intval('3.14159_') )
%
%is executed, the value v.x is the interval with left bound 3.14158 and
%right bound 3.14160 (correctly rounded) and gradient value v.dx = [0 0]
%(v is treated like a constant). Similarly,
%
%  u(2) = -4711
%
%produces u.x = [ -3 ; -4711 ] and u.dx = [ 1 0 ; 0 0 ]
%
%A simple two-dimensional Newton iteration for the Rosenbrock function
%   f = inline(' [ 400*x(1)*(x(1)*x(1)-x(2)) + 2*(x(1)-1) ; 200*x(1)*(x(1)*x(1)-x(2)) ]')
%starting at x = [1.1;0.5] is
%  x = gradientinit([1.1;0.5]);  y = f(x); x = x - y.dx\y.x
%producing
%gradient value x.x = 
%    1.0000
%    0.9255
%gradient derivative(s) x.dx = 
%     1     0
%     0     1
%
%For other examples, see demogradient and intlablogo.
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump
% modified 09/16/00     S.M. Rump  typo
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    handling of derivatives of sparse arrays
%                                    sparsegradient added
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/03/12     S.M. Rump  GradientSparseArrayDeriv and see removed
%

  if nargin==0
    setappdata(0,'INTLAB_GRADIENT_NUMVAR',0);
    return
  end
  
  if ~isa(v,'double') & ~isa(v,'intval')
    error('invalid initialization of dependent gradient variables')
  end
  setappdata(0,'INTLAB_GRADIENT_NUMVAR',prod(size(v)));
  dummy.init = v;
  r = gradient( dummy , 'gradientinit' );
