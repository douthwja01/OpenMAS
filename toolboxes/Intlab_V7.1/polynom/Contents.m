%INTLAB polynom toolbox
%
%Polynomial constructors
%  polynom      - polynomial constructor
%  intval       - Change to intval coefficients
%  subsasgn     - Subscripted assignment A(i) = 1
%  subsref      - Subscripted reference r = A(4) or p{3.1}
%
%Display of intervals (rigorous)
%  display      - Command window display of polynomial (use user defined default)
%  disp         - Display function for pop-up windows in debugger
%  disp_        - Display coefficients with uncertainty
%  infsup       - Display coefficients infimum and supremum
%  midrad       - Display coefficients midpoint and radius
%
%Interval arithmetic operations
%  plus         - Plus                              +
%  uplus        - Unary plus                        +
%  minus        - Minus                             -
%  uminus       - Unary minus                       -
%  mtimes       - Polynomial multiplication         *
%  mrdivide     - Polynomial right division         /
%  mldivide     - Polynomial left division          \
%  rdivide      - Polynomial right division         ./
%  ldivide      - Polynomial left division          .\
%  mpower       - Polynomial power                  ^
%
%Other polynomial operations
%  ctranspose   - Polynomial first derivative       '
%  pderiv       - Polynomial (higher) derivatives
%  deconv       - Polynomial division with remainder
%  abs          - Polynomial of interval absolute value coefficients
%  mag          - Polynomial of (real) absolute value coefficients
%  inf          - Polynomial of infimum of coefficients
%  sup          - Polynomial of supremum of coefficients
%  mid          - Polynomial of midpoint of coefficients
%  rad          - Polynomial of radius of coefficients
%  diam         - Polynomial of diameter of coefficients
%  real         - Polynomial of real part of coefficients
%  imag         - Polynomial of imaginary part of coefficients
%  qdist        - Metric distance
%  relerr       - Relative error of coefficients of two polynomials
%
%Utility routines
%  isreal       - Polynomial coefficients are real
%  isnan        - Polynomial coefficients contain NaN
%  isintval     - Polynomial coefficients are intervals
%  isfinite     - Polynomial coefficients are finite
%  isinf        - Logical result: true, if at least one coefficient not finite
%  vector       - Vector of coefficients of univariate polynomial 
%  find         - Find indices of nonzero coefficients
%  nnz          - Number of nonzero coefficients
%  end          - Last index: not supported because of ambiguouity
%
%Structural operations
%  degree       - Degree of polynomial
%  numvars      - Number of variables 
%  P(i)         - i-th coefficient of polynomial
%  P{x}         - Polynomial evaluation
%  polyval      - Polynomial evaluation
%  pshift       - Shift polynomial by constant
%  pscale       - Scale polynomial coefficients
%  ptrans       - Affine transformation of polynomial
%  removevars   - remove superfluous variables
%
%Auxiliary routines
%
%  roots        - Approximate roots of polynomial
%  rootbound    - Bound for all roots
%  rootsep      - root separation
%  plotpoly     - Plot of (interval) polynomial in <=2 variables
%  sylvester    - Sylvester matrix of two polynomials
%  permvars     - Permute variables of multivariate polynomials
%  bernsteincoeff - Bernstein coefficients
%  
%Generation of polynomials
%
%  randpoly     - univariate or multivariate random polynomial
%  bernstein    - Bernstein polynomials
%  chebyshev    - Chebyshev polynomial
%  chebyshev2   - Chebyshev polynomial of second kind
%  gegenbauer   - Gegenbauer polynomial
%  hermite      - Hermite polynomial
%  jacobi       - Jacobi polynomial
%  laguerre     - Laguerre polynomial
%  legendre     - Legendre polynomial
%  bugeaud      - Polynomial with very close roots
%
%Polynomial Comparison
%  eq           - Equal                             ==
%  ne           - Not equal                         ~=
%
%Verification routines and auxiliary
%  verifypoly   - Inclusion of simple and clustered roots of polynomial
%  typeadj      - Type adjustment
%  typeof       - Type for type adjustment
%
%Initialization of INTLAB intval package and system variables
%  polynominit  - Initialization and definition of defaults
%
%Demonstration, samples
%  demopolynom  - Some examples for using INTLAB polynom package
%

% written  10/04/02     S.M. Rump  Version 4

% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
