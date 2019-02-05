%INTLAB interval slope toolbox
%
%Slope constructors
%  slopeinit    - Initialization of dependent variables
%  horzcat      - Horizontal concatenation          [ , ]
%  vertcat      - Vertical concatenation            [ ; ]
%  subsasgn     - Subscripted assignment A(i,:) = 1
%  subsref      - Subscripted reference r = A(3,4)
%  slope        - Double to slope
%
%Display of slopes (rigorous)
%  display      - Command window display of slope
%  disp         - Display function for pop-up windows in debugger
%  infsup       - Display infimum and supremum
%  midrad       - Display midpoint and radius
%  disp_        - Display in "_" notation
%  slopeplot    - Plot of one-dimensional function and slope expansion
%
%slope arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Matrix multiply                   *
%  times        - Elementwise multiply              .*
%  mrdivide     - Slash or right division           /
%  rdivide      - Elementwise right division        ./
%  mpower       - Matrix power                      ^
%  power        - Elementwise power                 .^
%
%Other slope operations
%  transpose    - conjugate transpose '
%  ctranspose   - transpose .'
%  abs          - Absolute value
%  trace        - Trace
%  sum          - Sum
%  prod         - Product
%  newton       - Data for Newton iteration
%
%Utility routines
%  isnan        - True for Not a Number
%  isreal       - slope is real (for completeness)
%  isintval     - slope is interval (for completeness)
%  isfinite     - Interval is finite
%  isinf        - Interval is infinite
%  isempty      - Slope is empty in Matlab sense, i.e. []
%  issparse     - Slope has sparse structure
%  find         - find indices of nonzero elements
%  end          - determine last index
%
%Structural operations
%  band         - Extract band
%  diag         - Extract diagonal
%  tril         - Extract lower trianglar
%  triu         - Extract upper trianglar
%  length       - Length
%  size         - Size
%  dim          - Dimension of square matrix
%  sparse       - Type cast to sparse slope
%  full         - Type cast to full slope
%  reshape      - Reshape
%  repmat       - Duplicate arrays
%
%slope trigonometric functions
%  sin          - Sine
%  cos          - Cosine
%  tan          - Tangent
%  cot          - Cotangent
%  sec          - Secant
%  csc          - Cosecant
%  asin         - Inverse sine
%  acos         - Inverse cosine
%  atan         - Inverse tangent
%  acot         - Inverse cotangent
%  asec         - Inverse secant
%  acsc         - Inverse cosecant
%  sinh         - Hyperbolic sine
%  cosh         - Hyperbolic cosine
%  tanh         - Hyperbolic tangent
%  coth         - Hyperbolic cotangent
%  asinh        - Inverse hyperbolic sine
%  acosh        - Inverse hyperbolic cosine
%  atanh        - Inverse hyperbolic tangent
%  acoth        - Inverse hyperbolic cotangent
%
%slope exponential functions
%  exp          - Exponential
%  log          - Natural logarithm
%  log10        - Logarithm to base 10
%  sqr          - Square
%  sqrt         - Square root
%
%slope comparison of ".c" part
%  eq           - Equal                             ==
%  ne           - Not equal                         ~=
%  gt           - Greater than                      >
%  ge           - Greater than or equal             >=
%  lt           - Less than                         <
%  le           - Less than or equal                <=
%
%Verification routines and auxiliary
%  typeadj      - Type adjustment
%  typeof       - Type for type adjustment
%
%Initialization of INTLAB slope package and system variables
%  slopeinit    - Initialization and definition of INTLAB switches, also
%                   initialization of dependent variables
%
%Demonstration, samples
%  demoslope    - Some examples for using INTLAB slope package
%  demotest     - A test function for demoslope
%
%
%The slope package uses standard slope techniques (see Neumaier's book) with
%improvements for convex/concave functions and successive one-dimensional
%slopes according to
%  S.M. Rump: Expansion and Estimation of the Range of Nonlinear Functions,
%    Math. Comp. 65(216), pp. 1503-1512 (1996).
%

% written  12/06/98     S.M. Rump
%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
