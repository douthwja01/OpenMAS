%INTLAB some utility routines
%
%Matrix routines
%  isspd        - symmetric (Hermitian) is spd
%  isregular    - proof regularity for extremely ill-conditioned matrices
%  compmat      - Ostrowski's comparison matrix
%  gregk...     - Some sparse test matrices from Gregory/Karney test set
%  circulant    - Circulant out of first row
%  boothroyd    - very ill-conditioned matrix with known inverse
%
%Random numbers and matrices
%  randint      - Random integers in specified range
%  random       - Random numbers uniformly distrubuted in [-1,1] or
%                   within specified range
%  randsym      - Random symmetric matrix
%  randherm     - Random Hermitian matrix
%  randomc      - Complex random numbers with real and imaginary part as random
%  randmat      - random matrix with specified condition number
%  randorth     - random orthogonal matrix
%  randsvd      - random matrix with geometrically distributed singular values
%
%Floating-point related
%  pred         - Predecessor
%  succ         - Successor
%  subrealmin   - Smallest unnormalized positive floating-point number
%
%Other routines
%  format       - change interval output format
%  finish       - makes sure rounding mode is set to nearest after exiting Matlab
%  helpp        - intelligent help
%  binom        - Binomial coefficient (vectowise)
%  relerr       - Relative error
%  perms_       - Fast version of Matlab perms
%  sqr          - Square
%  odd          - logical integer odd
%  even         - logical integer even
%  factors      - List of factors of an integer
%  distbits     - Distance in bits
%  fletcher     - Sample routine for nonlinear systems
%

% written  11/30/98     S.M. Rump
% modified 12/15/01     S.M. Rump  Routine Fletcher added
% modified 11/23/05     S.M. Rump  band and bandwidth removed (already in /intval)
%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
