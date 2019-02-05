%Reference implementations for accurate summation and dot product algorithms
%
%Generation of extremely ill-conditioned sums and dot products
%  GenSum       - generation of extremely ill-conditioned sums
%  GenDot       - generation of extremely ill-conditioned dot products
%
%New summation, dot product and related routines: high accuracy
%  FastAccSum    - faithful rounding of sum(p_i)
%  AccSum        - faithful rounding of sum(p_i)
%  AccSumHugeN   - faithful rounding of sum(p_i) for large dimension
%  AccSumK       - K-fold faithful rounding of sum(p_i)
%  AccSign       - sign of sum(p_i)
%  PrecSum       - accurate and fast up to large condition number
%  FastPrecSum   - accurate and fast up to large condition number
%  NearSum       - sum(p_i) rounded to nearest
%  DownSum       - sum(p_i) rounded downwards
%  UpSum         - sum(p_i) rounded upwards
%
%
%New summation, dot product and related routines: high precision
%  Sum2          - summation with quad precision
%  SumK          - summation with K-fold precision
%  SumKL         - dot product computed in K-fold precision, result stored in L parts
%  Sum_          - summation in K-fold precision
%  Dot2          - dot product with quad precision
%  Dot2Err       - dot2 with rigorous error bounds without directed rounding
%  DotK          - dot product with K-fold precision
%  DotKL         - dot product with K-fold precision, result stored in L parts
%  Dot_          - easy to use dot product for vectors and matrices
%
%
%Error-free transformations
%  TwoSum        - transformation of a+b into x+y with x=fl(a+b)
%  FastTwoSum    - as TwoSum provided input is ordered in absolute value
%  TwoProduct    - transformation of a*b into x+y with x=fl(a*b)
%  Split         - transformation of a into two 'halves' x+y
%  ExtractVector - extract higher and lower part of vector
%  Transform     - Transformation of vector into high part and low order vector
%
%
%Utility routines for new summation and dot product routines
%  NextPowerTwo  - next power of 2 of integer without branch
%  ufp           - unit in the first place
%
%
%Reference implementations of competitors
%  SumXBLAS      - summation as in XBLAS
%  DotXBLAS      - dot product as in XBLAS
%  PriestSum     - Priest's doubly compensated summation
%  ZDSum         - Zielke/Drygalla summation
%
%
%Application program
%  InvIllco      - Inverse of extremely ill-conditioned matrices
%

% written  10/27/08     S.M. Rump
%

%
%New algorithms based on
%
%Accurate summation and dot product with specified precision
%  T. Ogita, S.M. Rump, and S. Oishi. Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005.
%
%Accurate summation and dot product with specified accuracy
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, SIAM J. Sci. Comput., 31(1):189-224, 2008.
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%
%Fast and high precision summation and dot products
%  S.M. Rump, T. Ogita, and S. Oishi. Fast high precision summation. 
%    Nonlinear Theory and Its Applications (NOLTA), IEICE, 1(1), 2010. 
%    [received the "NOLTA Best Paper Award" by the IEICE Engineering Sciences
%    Society]. 
%

%
% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
