function r = taylorinit(v,k)
%TAYLORINIT   Initialization of taylor variable
%
%  x = taylorinit(v,k)
%
%The dependent variable x is identified and initialized with v, and Taylor
%coefficients up to k>0 are calculated. For example,
%
%  u = taylorinit(-3,5)
%
%initializes "u" to be a Taylor variable with value -3 and Taylor coefficients
%up to order k=5, i.e. [-3,1,0,0,0,0]. The default for k is 4. After initialization,
%the order of Taylor coefficients can be obtained by
%
%  k = taylororder
%
%In our example, the result will be k=5. The evaluation
%
%  f = inline('x*sin(x)-exp(x^2)')
%  y = f(u)
%
%computes the Taylor coefficients of f up to order 5. Note the first component of y
%is f(u), followed by f'(u), f''(u)/2!, f'''(u)/3!, etc. For example,
%
%  t0 = y{0}
%
%is the 0-th Taylor coefficients, namely f(u), and 
%
%  sum( y{0:k} .* e.^(0:k) )
%
%is approximately equal to f(u+e) up to O(e^(k+1)). The subsequent statement
%
%  v = taylor( intval('3.14159_') )
%
%initializes v to be the constant with interval value infsup(3.14158,3.14160).
%
%Taylor expansions are implemented for univariate functions. However, it is 
%sometimes useful to evaluate a function at an array of points. The statement
%
%  w = taylorinit((-3:1)',5)
%
%initializes w to be a vector of dependent variables with values -3,-2,-1,0,1. All 
%operations with u are executed pointwise, so y=f(w) is executed as w.*sin(w)-exp(w.^2).
%In this case
%
%  y(3)
%
%is the Taylor expansion of f at -1 up to order 5.
%
%For other examples, see demotaylor.
%

% written  05/21/09     S.M. Rump
% modified 09/15/10     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  setappdata(0,'INTLAB_TAYLOR_END',0);

  if nargin==0
    setappdata(0,'INTLAB_TAYLOR_ORDER',0);
    return
  end
  
  if ~isa(v,'double') & ~isa(v,'intval')
    error('invalid initialization of dependent Taylor variables')
  end
  if nargin==1
    setappdata(0,'INTLAB_TAYLOR_ORDER',4);
  else
    setappdata(0,'INTLAB_TAYLOR_ORDER',k);
  end
  dummy.init = v;
  r = taylor( dummy , 'taylorinit' );
