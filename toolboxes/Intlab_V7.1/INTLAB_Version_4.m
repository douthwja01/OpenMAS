%Help file for INTLAB Version 4
%
%Polynom toolbox added (for univariate and multivariate polynomials)
% $$$ Disclaimer: Sometimes slow due to interpretation overhead, expecially for interval polynomials $$$
%including basic operations, various utility programs (Chebyshev, Bernstein ...) and a verification routine 
%for multiple roots.
%
%Functionality of abs and abss had to be interchanged. The natural extension of abs for
%  interval argument seems to be the function yielding the **interval of all absolute values**,
%  whereas the maximum absolute value of numbers within the interval, a real number, is now
%  computed by abss. Be sure to change existing programs according to that.
%
%Additions/changes:
%
%  - improved residual approximation and inclusion "lssresidual" added with corresponding
%      switch in intvalinit; implies improved inclusion for linear systems
%  - randmat added: random matrix with specified condition number
%  - displaywidth added to adjust output (see intvalinit)
%  - spdiags for interval matrices added
%  - plotintval added (try it, may give nice pictures)
%  - some utility functions added such as isinf, isfinite, ...
%  - some functions for sparse structures and others added such as nnz, nonzeros, find, ...
%  - functionality of functions sparse.m improved
%  - interval power redesigned, especially for integer exponents including point intervals being integer
%  - verifylss etc. deliver intval(NaN) in case of non-success
%  - setround to work under several OS without change
%  - to avoid confusion, all INTLAB global variables now begin with INTLAB_
%  - functions isinf and isfinite added
%  - index "end" now possible for INTLAB data types
%  - some changes due to different behavior of Matlab 6.5
%
%  Others:
%
%  - some performance improvements, especially for sparse data
%  - some other utility functions added
%
