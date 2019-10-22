%INTLAB interval hessian toolbox
%
%Hessian constructors
%  hessianinit  - Initialization of dependent variables
%  horzcat      - Horizontal concatenation          [ , ]
%  vertcat      - Vertical concatenation            [ ; ]
%  subsasgn     - Subscripted assignment A(i,:) = 1
%  subsref      - Subscripted reference r = A(3,4)
%  hessian      - Double to gradient
%
%Display of hessians and interval hessians (rigorous)
%  display      - Command window display of hessian
%  disp         - Display function for pop-up windows in debugger
%  realimag     - Real and Imaginary part separately
%  infsup       - Display infimum and supremum
%  midrad       - Display midpoint and radius
%  disp_        - Display in "_" notation
%
%Hessian arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Matrix multiply                   *
%  times        - Elementwise multiply              .*
%  mrdivide     - Slash or right division           /
%  mldivide     - Backslash or left division        \
%  rdivide      - Elementwise right division        ./
%  ldivide      - Elementwise left division         .\
%  mpower       - Matrix power                      ^
%  power        - Elementwise power                 .^
%  intersect    - Intersection
%
%Other hessian operations
%  abs          - Absolute value
%  inf          - Infimum
%  sup          - Supremum
%  mid          - Midpoint
%  rad          - Radius
%  diam         - Diameter
%  real         - Real part
%  imag         - Imaginary part
%  trace        - Trace
%  sum          - Sum
%  prod         - Product
%  ctranspose   - Complex conjugate transpose       '
%  transpose    - Transpose                         .'
%
%Utility routines
%  isnan        - True for Not a Number
%  isreal       - Hessian is real
%  isintval     - Hessian is intval
%  isfinite     - Interval is finite
%  isinf        - Interval is infinite
%  isempty      - Hessian is empty in Matlab sense, i.e. []
%  emptyintersect - detect empty intersections
%  issparse     - Hessian has sparse structure
%  find         - find indices of nonzero elements
%  all          - Determine if all array elements are nonzero
%  any          - Determine if any array elements are nonzero
%  logical      - Convert hessian values to logical
%  nnz          - number of nonzero elements
%  end          - determine last index
%
%Structural operations
%  full         - Convert to full hessian
%  sparse       - Convert to sparse hessian
%  band         - Extract band
%  diag         - Extract diagonal
%  tril         - Extract lower triangular
%  triu         - Extract upper triangular
%  bandwidth    - Bandwidth
%  length       - Length
%  size         - Size
%  dim          - Dimension of square matrix
%  reshape      - Reshape
%  repmat       - Duplicate arrays
%
%Hessian trigonometric functions (rigorous, real and complex)
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
%Hessian exponential functions (rigorous, real and complex)
%  exp          - Exponential
%  log          - Natural logarithm
%  log10        - Logarithm to base 10
%  sqr          - Square
%  sqrt         - Square root
%
%Hessian comparison of ".x" part
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
%Initialization of INTLAB hessian package and system variables
%  hessianinit  - Initialization and definition of INTLAB switches, also
%                   initialization of dependent variables
%
%Demonstration, samples
%  demohessian - Some examples for using INTLAB hessian package
%
%
%The hessian package uses forward mode of automatic differentiation. For
%an introduction to forward differentation, see
%  L.B. Rall: Automatic Differentiation: Techniques and Applications,
%    Lecture Notes in Computer Science 120, Springer, 1981.
%

% written  04/04/04     S.M. Rump  INTLAB Version 5
%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
