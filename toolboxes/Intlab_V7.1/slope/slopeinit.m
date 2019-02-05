function x = slopeinit(xs,X)
%SLOPEINIT    Initialization of slope expansion
%
%  x = slopeinit(xs,X)
%
%This statement with interval quantities xs and X followed by
%  a function evaluation
%
%  y = f(x)
%
%using arithmetic operations and standard functions (overloaded by
%slope operators) implies 3 assertions:
%
%  1)  f(xs) in y.c
%  2)  f(X)  in y.r
%  3)  For all xp in X there exists s in y.s with  f(xp) = f(xs) + s*(xp-xs) .
%
%The expansion point xs and expansion interval X must be real quantities,
%  non-interval input data is converted into interval data to ensure
%  correctness of 1..3).
%
%The call
%
%  slopeinit
%
%without input parameters initializes the slope expansion to xs = X = [].
%
%The slope of arrays is stored in the 'next dimension'. So y.s is 
%  3-dimensional for a slope row vector y. Since Matlab does not 
%  support multi-dimensional sparse arrays, y.s is not accessible 
%  if y has more than one column.
%
%As a simple example of slopes
%
%  xs = [ -3 ; 7.15 ] ;  X = intval(' [ -3.1,-2.9 ]  [ 7,7.1 ] ') ;
%  u = slopeinit( xs , X );
%
%initializes the slope package to have two dependent variables. The
%expansion point is xs=[-3;7.15] with expansion interval X.
%
%For a function f from R^2 to R, y=f(u) computed by overloaded slope
%operators has access y.c, y.r and y.s with the properties listed above.
%Internally, slopes are computed and stored according to
%
%  S.M. Rump: Expansion and Estimation of the Range of Nonlinear Functions,
%    Math. Comp. 65(216), p. 1503-1512 (1996).
%
%In contrast to the usual definition of slopes (see Neumaier's book),
%
%  - rather than one n-dimensional slope, n one-dimensional slopes
%      are calculated treating the function
%        f(X_1,...,X_k-1,z,xs_k+1,...,xs_n)
%      as a function in one variable z.
%  - intersections for the two definitions of multiplication and division
%      sharpens results (note this is not possible for n-dimensional slopes)
%  - intersections of the computed range y.r with y.c+y.s(X-xs) sharpens
%      results
%  - the sharp formulas for slopes of convex and concave functions are
%      used (for details, see the paper).
%
%Slope are useful for verified inclusion of clustered or multiple zeros.
%Sometimes, the range estimation is better than naive interval arithmetic
%or gradient expansion.
%
%For a simple plot of one-dimensional slopes, see slopeplot.
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/06/04     S.M. Rump  handling of derivatives of sparse arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/03/12     S.M. Rump  SlopeSparseArrayDeriv removed
%

  if nargin==0
    INTLAB_SLOPE.NUMVAR = 0;
    INTLAB_SLOPE.Xxs = [];
    setappdata(0,'INTLAB_SLOPE',INTLAB_SLOPE);
    return
  end

  if ~isreal(xs) | ~isreal(X)
    error('Complex numbers not allowed in slope expansion')
  end

  if ( ~isa(xs,'double') & ~isa(xs,'intval') ) | ...
     ( ~isa(X,'double') & ~isa(X,'intval') )
    error('invalid initialization of slopes: input must be double or intval')
  end

  if ~isequal(size(xs),size(X))
    error('invalid initialization of slopes: dimensions do not match')
  end

  % initialize INTLAB constants
  INTLAB_SLOPE.NUMVAR = prod(size(xs));
  dummy.xs = intval(xs);
  dummy.X = intval(X);
  INTLAB_SLOPE.Xxs = dummy.X - dummy.xs;
  INTLAB_SLOPE.Xxs = INTLAB_SLOPE.Xxs(:);
  setappdata(0,'INTLAB_SLOPE',INTLAB_SLOPE);
  x = slope( dummy, 'slopeinit' );
