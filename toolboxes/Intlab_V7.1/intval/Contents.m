%INTLAB interval toolbox
%
%Interval constructors
%  intval       - Intval constructor
%  cintval      - Complex intval type cast
%  infsup       - Infimum/supremum to interval
%  midrad       - Midpoint/radius to interval
%  gradient     - Gradient to interval
%  slope        - Slope to interval
%  spones       - Sparse interval of ones
%  horzcat      - Horizontal concatenation          [ , ]
%  vertcat      - Vertical concatenation            [ ; ]
%  subsasgn     - Subscripted assignment A(i,:) = 1
%  subsref      - Subscripted reference r = A(3,4)
%  tocmplx      - replaced by cintval
%
%Display of intervals (rigorous)
%  display      - Command window display of interval (use user defined default)
%  disp         - Display function for pop-up windows in debugger
%  displaywidth - Set width of display
%  disp_        - Display with uncertainty
%  infsup       - Display infimum and supremum
%  midrad       - Display midpoint and radius
%  displayinner - Display of inner inclusion
%  realimag     - Real and Imaginary part separately
%  str2intval   - Rigorous conversion string to intval
%  disp2str     - Special display into common exponent and mantissa
%  plotintval   - Plots real or complex intervals
%
%Interval arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Matrix multiply                   *
%  times        - Elementwise multiply              .*
%  mldivide     - Backslash, linear system solving  \
%  mrdivide     - Slash or right division           /
%  rdivide      - Elementwise right division        ./
%  ldivide      - Elementwise left division         .\
%  mpower       - Matrix power                      ^
%  power        - Elementwise power                 .^
%  intersect    - Intersection
%  hull         - Hull or convex union
%
%Other interval operations
%  setround     - Set rounding mode
%  getround     - Get rounding mode
%  ctranspose   - Complex conjugate transpose       '
%  transpose    - Transpose                         .'
%  conj         - Conjugate
%  abs          - Interval absolute value, result real interval
%  mag          - Absolute value, result double
%  mig          - Mignitude, result double
%  inf          - Infimum
%  inf_         - Infimum (for problems with inf)
%  sup          - Supremum
%  mid          - Midpoint
%  rad          - Radius
%  diam         - Diameter
%  real         - Real part
%  imag         - Imaginary part
%  trace        - Trace
%  sum          - Sum
%  prod         - Product
%  norm         - Norm
%  qdist        - Metric distance
%  relerr       - Relative error of intervals
%  relacc       - Relative accuracy of intervals
%  gershgorin   - Gershgorin discs
%
%Utility routines
%  isnan        - True for Not a Number
%  isreal       - Interval is real
%  isintval     - for completeness
%  isfinite     - Interval is finite
%  isinf        - Interval is infinite
%  isempty      - Interval is empty in Matlab sense, i.e. []
%  emptyintersect  - Check for empty intersection
%  iszero       - Interval is zero (componentwise)
%  issparse     - Interval has sparse structure
%  spy          - Spy sparse interval matrix
%  find         - find indices of nonzero elements
%  all          - Determine if all array elements are nonzero
%  any          - Determine if any array elements are nonzero
%  logical      - Convert intval values to logical
%  nnz          - number of nonzero elements
%  nonzeros     - vector of nonzero elements
%  end          - determine last index
%
%Structural operations
%  band         - Extract band
%  diag         - diagonal
%  tril         - Extract lower triangular
%  triu         - Extract upper triangular
%  bandwidth    - Bandwidth
%  length       - Length
%  size         - Size
%  dim          - Dimension of square matrix
%  reshape      - Reshape
%  repmat       - Duplicate arrays
%  permute      - Permute array dimensions
%  ipermute     - Inverse permute array dimensions
%  sparse       - Type cast to sparse interval matrix
%  full         - Type cast to full interval matrix
%  spdiags      - Generate sparse matrix out of diagonals
%  toeplitz     - create interval toeplitz matrix
%  circulant    - create interval circulant matrix
%
%Interval trigonometric functions (rigorous, real and complex)
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
%  modpi2       - Precise mod pi/2 (for argument reduction)
%  stdfctsdata  - To generate data for rigorous standard functions
%
%Interval exponential functions (rigorous, real and complex)
%  exp          - Exponential
%  log          - Natural logarithm
%  log10        - Logarithm to base 10
%  sqr          - Square
%  sqrt         - Square root
%
%Interval higher transcendental functions (rigorous, real)
%  erf          - error function
%  erfc         - complementary error function
%
%Interval Comparison
%  eq           - Equal                             ==
%  ne           - Not equal                         ~=
%  gt           - Greater than                      >
%  ge           - Greater than or equal             >=
%  lt           - Less than                         <
%  le           - Less than or equal                <=
%  in           - Contained in
%  in0          - Contained in interior
%  iszero       - Interval is zero (componentwise)
%
%Verification routines and auxiliary
%  verifylss    - Verified linear system solver including
%                    rectangular and sparse systems
%  structure    - Specification of user-defined structured matrices
%  structlss    - Verified solution of structured linear systems
%  lssresidual  - Improved approximation and inclusion of residual
%  plotlinsol   - Solution of 2x2 interval linear systems
%  verifyeig    - Verified eigenvalue inclusion (simple and clusters)
%                    together with basis of invariant subspaces,
%                    ordinary and generalized eigenvalue problem
%  structeig    - Verified eigenvalue/vector inclusion (simple and clusters)
%                    for structured input matrix
%  verifyquad   - Verified quadrature
%  verifynlss   - Verified nonlinear system solver
%  verifynlss2  - Verified nonlinear system solver for double roots
%  test         - Sample nonlinear test functions
%  compmat      - Ostrowski's comparison matrix
%  inv          - Verified inverse of square matrix
%  typeadj      - Type adjustment
%  typeof       - Type for type adjustment
%
%Initialization of INTLAB intval package and system variables
%  intvalinit   - Initialization and definition of defaults
%  stdfctsinit  - Initialization, defaults for rigorous standard functions
%
%Demonstration, samples
%  demointval   - Some examples for using INTLAB intval package
%
%
%For introduction to interval arithmetic see, for example,
%  A. Neumaier: Interval Methods for Systems of Equations, Encyclopedia
%    of Mathematics and its Applications, Cambridge University press, 1990.
%
%INTLAB implementation of interval arithmetic is based on
%  S.M. Rump: Fast and Parallel Interval Arithmetic,
%    BIT 39(3), 539-560, 1999.
%and fast conversion of inf-sup to mid-rad based on
%  S. Oishi, S.M. Rump: Fast verification of solutions of matrix equations, 
%    Numerische Mathametik 90, 755-773, 2002.
%
%Real interval standard functions based on
%  S.M. Rump: Rigorous and portable standard functions,
%    BIT 41(3), 540-562, 2001.
%
%Complex interval standard functions based on
%  N.C. Boersken: Komplexe Kreis-Standardfunktionen, Freiburger
%    Intervallberichte 78/2.
%

% written  11/30/98     S.M. Rump
% modified 03/06/99     S.M. Rump  Version 2
% modified 11/11/99     S.M. Rump  Version 3
% modified 03/11/02     S.M. Rump  Version 3.1
% modified 12/04/05     S.M. Rump  Version 5.2
% modified 05/26/05     S.M. Rump  Version 5.3
% modified 05/31/05     S.M. Rump  Version 5.4
% modified 11/11/07     S.M. Rump  Version 5.5
% modified 03/30/10     S.M. Rump  Version 6
% modified 10/16/12     S.M. Rump  Version 7
% modified 06/24/13     S.M. Rump  Version 7.1
%

% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
