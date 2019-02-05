% INTLAB  INTerval LABoratory. 
% Version 7.1   June-24-2013   966 files, 849 m-files, 32,346 lines of Matlab code w/o comments 
%                 (54,324 with comments),
%             supported by a test-case library with more than 60,000 lines of Matlab code w/o comments
%
%===> For license information and terms of use, see README.TXT.
%
%
%Directories and toolboxes
%  intval       - Interval package
%  gradient     - Automatic differentiation package
%  hessian      - Automatic Hessian package
%  taylor       - Automatic Taylor package
%  slope        - Automatic slope package
%  polynom      - Polynomial package (univariate and multivariate)
%  long         - Rudimentary long package
%  utility      - Some useful functions
%  AccSum       - Reference implementations for sum/dot routines
%  demos        - several demo routines
%
%Auxiliary files
%  intlablogo   - INTLAB logo display
%  startintlab  - Initialization of INTLAB
%  startup      - Calls startintlab: adapt to your local installation
%  intlabsetting   - Current setting of INTLAB control variables
%  INTLAB_Version_2     - Additions and changes in version 2
%  INTLAB_Version_3     - Additions and changes in version 3
%  INTLAB_Version_3.1   - Additions and changes in version 3.1
%  INTLAB_Version_4     - Additions and changes in version 4
%  INTLAB_Version_4.1   - Additions and changes in version 4.1
%  INTLAB_Version_4.1.1 - Additions and changes in version 4.1.1
%  INTLAB_Version_4.1.2 - Additions and changes in version 4.1.2
%  INTLAB_Version_5     - Additions and changes in version 5
%  INTLAB_Version_5.1   - Additions and changes in version 5.1
%  INTLAB_Version_5.2   - Additions and changes in version 5.2
%  INTLAB_Version_5.3   - Additions and changes in version 5.3
%  INTLAB_Version_5.4   - Additions and changes in version 5.4
%  INTLAB_Version_5.5   - Additions and changes in version 5.5
%  INTLAB_Version_6     - Additions and changes in version 6
%  INTLAB_Version_7     - Additions and changes in version 7
%  INTLAB_Version_7.1   - Additions and changes in version 7.1
%  FAQ          - Frequently asked questions
%  Readme       - Installation, a little tutorial and miscellaneous
%
%
%INTLAB based on
%  S.M. Rump: INTLAB - INTerval LABoratory, in "Developments in
%    Reliable Computing", ed. T. Csendes, Kluwer Academic Publishers,
%    77-104, 1999.
%
%INTLAB implementation of interval arithmetic is based on
%  S.M. Rump: Fast and Parallel Interval Arithmetic,
%    BIT 39(3), 539-560, 1999.
%and
%  S. Oishi, S.M. Rump: Fast verification of solutions of matrix equations, 
%    Numerische Mathametik 90, 755-773, 2002.
%
%Real interval standard functions based on
%  S.M. Rump: Rigorous and portable standard functions,
%    BIT 41(3), 540-56, 2001.
%
%Complex interval standard functions based on
%  N.C. Boersken: Komplexe Kreis-Standardfunktionen, Freiburger
%    Intervallberichte 78/2.
%
%Accurate summation and dot product with specified _precision_
%  T. Ogita, S.M. Rump, and S. Oishi. Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005.
%  S.M. Rump, T. Ogita, and S. Oishi. Fast High Precision Summation. 
%     Nonlinear Theory and Its Applications (NOLTA), IEICE, 1(1), 2010.
%  S.M. Rump: Ultimately Fast Accurate Summation, SIAM Journal on 
%     Scientific Computing (SISC), 31(5):3466-3502, 2009.
%
%Accurate summation and dot product with specified _accuracy_
%  S.M. Rump, T. Ogita, and S. Oishi: Accurate Floating-point Summation I: 
%     Faithful Rounding. SIAM Journal on Scientific Computing (SISC), 
%     31(1): 189-224, 2008.
%  S.M. Rump, T. Ogita, and S. Oishi: Accurate Floating-point Summation II: 
%     Sign, K-fold Faithful and Rounding to Nearest. SIAM Journal on Scientific 
%     Computing (SISC), 31(2):1269-1302, 2008.
%  S.M. Rump: Ultimately Fast Accurate Summation, SIAM Journal on 
%     Scientific Computing (SISC), 31(5):3466-3502, 2009
%
%The linear system solver including inner inclusions is completely redesigned and 
%now based on
%  S.M. Rump. Accurate solution of dense linear systems, Part II: Algorithms using 
%    directed rounding. Journal of Computational and Applied Mathematics (JCAM), 
%    242:185–212, 2013.
%
%Least squares problems and underdetermined linear systems are based on
%  S.M. Rump. Verified Bounds for Least Squares Problems and Underdetermined Linear Systems. 
%    SIAM J. Matrix Anal. Appl. (SIMAX), 33(1):130–148, 2012.
%and new componentwise error estimates on
%  S.M. Rump: Improved componentwise verified error bounds for least squares problems 
%    and underdetermined linear systems, to appear.
%
%For references to verification algorithms cf. the corresponding routines. Many
%of them are discussed in 
%  S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic.
%    Acta Numerica, 19:287-449, 2010. 
%
%Slopes based on
%  R. Krawzcyk, A. Neumaier: Interval slopes for rational functions
%    and associated centered forms, SIAM J. Numer. Anal. 22, 604-616
%    (1985)
%with improvements based on
%  S.M. Rump: Expansion and Estimation of the Range of Nonlinear Functions,
%    Math. Comp. 65(216), pp. 1503-1512, 1996.
%
%Gradients, Hessians and Taylor expansion based on standard forward mode, cf. for example
%  L.B. Rall: Automatic Differentiation: Techniques and Applications,
%    Lecture Notes in Computer Science 120, Springer, 1981.
%
%Long floating point and interval arithmetic based on standard techniques
%  with adaptation and speed up for interpretative systems.
%

% written  12/30/98     S.M. Rump
% modified 03/06/99     S.M. Rump  Version 2
% modified 11/16/99     S.M. Rump  Version 3
% modified 03/07/02     S.M. Rump  Version 3.1
% modified 12/08/02     S.M. Rump  Version 4
% modified 12/27/02     S.M. Rump  Version 4.1
% modified 01/22/03     S.M. Rump  Version 4.1.1
% modified 11/18/03     S.M. Rump  Version 4.1.2
% modified 04/04/04     S.M. Rump  Version 5
% modified 06/04/05     S.M. Rump  Version 5.1
% modified 12/20/05     S.M. Rump  Version 5.2
% modified 05/26/06     S.M. Rump  Version 5.3
% modified 05/31/07     S.M. Rump  Version 5.4
% modified 11/05/08     S.M. Rump  Version 5.5
% modified 05/08/09     S.M. Rump  Version 6
% modified 12/12/12     S.M. Rump  Version 7
% modified 06/24/13     S.M. Rump  Version 7.1
%

% Copyright (c) Siegfried M. Rump, head of the Institute for Reliable Computing, 
%               Hamburg University of Technology
%
% For license information and terms of use, see README.TXT.
% 
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% ::  Some references to papers using INTLAB are collected on the
% ::    INTLAB homepage   http://www.ti3.tuhh.de/ 
% ::  If you have additional references to add, please send me
% ::  a mail ( rump [at] tuhh.de ) 
%