%% DEMOINTVAL  Interval computations in INTLAB
%

%%
% A key to interval operations in INTLAB is changing the rounding
% mode. Following we ensure "rounding to nearest".

setround(0)                    

%% How to define an interval I
% There are four possibilities to generate an interval, the first is a simple
% typecast of a real or complex quantity, for example a matrix.
% It uses Matlab 
% conversion, i.e. the first component does *not* necessarily contain "2.3".
% This is because Matlab first converts "2.3" into binary format, and then
% the type cast is called.

format compact long infsup
u = intval( [ 2.3 -4e1 ; 3 0 ] )
      
%% How to define an interval II
% The second possibility is to use INTLAB conversion of constants. 
% In this case the argument is a string and INTLAB verified conversion
% to binary is called rather than Matlab conversion. 
      
u = intval( '0.1 -3.1 5e2 .3e1' )

%%
% The first component, for example, definitely contains "0.1". Since 0.1
% is not finitely representable in binary format, the radius of the first
% component must be nonzero.

u.rad
     
%%
% Generating an interval by an input string always produces a column vector.
% To change "u" into a 2 x 2 matrix, use

reshape(u,2,2)

%%
% Note that arrays are stored columnwise in Matlab.      
    
%% How to define an interval III
% The third possibility is by specification of midpoint and radius.

u = midrad( [ -3e1+2i ; .4711 ; 3 ] , 1e-4 )
        
%% How to define an interval IV
% The fourth possibility is by specification of infimum and supremum
% 

u = infsup( [ 1 2 ] , [ 4 7 ] ) 
  
%% Output formats of intervals I
% The output format can be changed using the different Matlab formats,
% for example
             
format midrad long e
X = midrad( [ -3e1+2i ; .4711 ; 3 ] , 1e-4 )

%%
% or

format infsup long
X

%%
% or with a common exponent

1e4*X

%% Output formats of intervals II
% With two arguments the functions infsup and midrad define an interval,
% with one argument they control the output of an interval:
  
format short
u = infsup( [ 1 2 ] , [ 4 7 ] ); 
infsup(u), midrad(u)
           
%% Rigorous output
% Note that output in INTLAB is rigorous. That means, 
% 
%   left <= ans <= right     for inf/sup notation
%   ans  in  mid+/-rad       for mid/rad notation
%
% where  ans  is the true (real or complex) answer, and left,right,
% mid,rad are the numbers corresponding to the _displayed_ decimal figures.


%% Output formats of intervals III
% A convenient way to display narrow intervals is the following:
      
x=midrad(pi,1e-8);
format short, infsup(x), midrad(x), disp_(x)
format long, infsup(x), midrad(x), disp_(x)      
format short          

%%
% Mathematically the following is true: Form an interval of the displayed midpoint and
% a radius of 1 unit of the last displayed decimal figure, then this is a correct inclusion
% of the stored interval.


%% Changing interval output permanently   
% The interval output format can be changed permanently, for example,
% to infimum/supremum notation: 
   
u = midrad( [ -3e1+2i ; .4711 ; 3 ] , 1e-4 );
format infsup
u
      
%%
% or to midpoint/radius notation: 

format midrad
u
      
%%
% or to display with uncertainties depicted by "_": 
  
format _
u


%% Display with uncertainty
% Display with uncertainty makes only sense for sufficiently narrow intervals.
% If the accuracy becomes too poor, INTLAB changes automatically to inf-sup
% or mid-rad display for real or complex intervals, respectively:

for k=-5:-1
  disp_(midrad(pi,10^k))
end

%% Newton iteration
% The following code is an interval Newton iteration to include sqrt(2).

format long
f = @(x)(x^2-2);                                % Function f(x) = x^2-2
X = infsup(1.4,1.7);                            % Starting interval
for i=1:4
  xs = X.mid;                                   % Midpoint of current interval
  Y = f(gradientinit(X));                       % Gradient evaluation of f(X)
  X = intersect( X , xs-f(intval(xs))/Y.dx )    % Interval Newton step with intersection
end

%%
% The "display_" output format shows nicely the quadratic convergence. 
%
% The last displayed result (which is in fact an interval) proves that the true
% value of sqrt(2) is between 1.41421356237308 and 1.41421356237310.
% Indeed, sqrt(2)=1.41421356237309504...
%
% The format "long e" in Matlab displays the most figures. With this we see that
% the internal accuracy of the final X is in fact even better, the width is
% only 2 units in the last place.

format long e
X
format short
          

%% Invoking interval operations
% An operation uses interval arithmetic if at least one of the operands is of type
% intval. For example, in 
   
u = intval(5); 
y = 3*u-exp(u)

%%
% the result y is an inclusion of 15-exp(5). However, in
       
u = intval(5); 
z = 3*pi*u-exp(u)

%%
% the first multiplication "3*pi" is a floating point multiplication. Thus
% it is not guaranteed that the result z is an inclusion of 15pi-exp(5).

%% Interval matrix operations
% INTLAB is designed to be fast. Case distinctions in interval multiplication
% can slow down computations significantly due to interpretation overhead. 
% Therefore, there is a choice between 'fast' and 'sharp' evaluation of interval
% matrix products. This applies only to 'thick' intervals, i.e. intervals with
% nonzero diameter. 
%

%% Sharp interval multiplication
% In the following example, c is a real random matrix, C is an interval matrix
% with diameter zero (a thin interval matrix), and CC is an interval matrix with
% nonzero diameter (a thick interval matrix), all of dimension nxn for n=1000. 
% First we measure the computing time with option 'SharpIVmult'.

n = 1000;
c = randn(n); 
C = intval(c);
C_ = midrad(c,.1);
intvalinit('SharpIVmult')
tic, scc = c*c; toc
tic, sCC = C*C; toc
tic, sCC = C*C_; toc
tic, sCC__ = C_*C_; toc

%%
% As can be seen, there is not much penalty if not both matrices are thick
% interval matrices; then, however, computation is slowed down significantly.

%% Fast interval multiplication

intvalinit('FastIVmult')
tic, fcc = c*c; toc
tic, fCC = C*C; toc
tic, fCC = C*C_; toc
tic, fCC__ = C_*C_; toc
max(max(diam(fCC__)./diam(sCC__)))

%%
% As can be seen there is again not much penalty if not both matrices are thick. 
% However, the 'fast' implementation is much faster than the 'sharp' at the cost of
% a little wider output. If intervals are very wide and any overestimation
% cannot be afforded (as in global optimization), the option 'SharpIVmult'
% is recommended. It is shown in
%
% S.M. Rump: Fast and parallel interval arithmetic. BIT Numerical Mathematics, 
% 39(3):539-560, 1999
%
% that the maximum (componentwise) overestimation by the option 'FastIVmult'
% compared to 'SharpIVmult' is a factor 1.5, for real and complex intervals.

%% Acceleration by vector/matrix notation
% It is advisable to use vector/matrix notation when using interval operation.
% Consider

n = 1000; x = 1:n; y = intval(x);
tic
for i=1:n
  y(i) = y(i)^2 - y(i);
end
t1 = toc

%%
% This simple initialization takes considerable computing time. Compare to

tic
y = intval(x);
y = y.^2 - y;
t2 = toc
ratio = t1/t2

%%
% Sometimes code looks more complicated, a comment may help. It is worth it.

%% Overestimation of interval operations
% Note that the main principle of interval arithmetic is that for given intervals
% A,B and an operation o, the result a o b is included in the interval result A o B
% for all a in A and all b in B. Since the result must be an interval, overestimations
% cannot be avoided in many situations. For example, in

close, kmax = 40; i = sqrt(-1); a=midrad(2,.7); b=midrad(1-i,1);
plotintval(3*a-exp(b)); hold on
phi = linspace(0,2*pi,kmax);
[A,B] = meshgrid( mid(a)+rad(a)*exp(i*phi) , mid(b)+rad(b)*exp(i*phi) );
plot(3*A-exp(B))
hold off

%%
% the blue circle is the result of the interval operations, whereas the many
% circles approximate the power set operation (see also the INTLAB demo). Another
% reason for overestimation are dependencies, see below. 

%% Interval standard functions               
% Interval standard functions in INTLAB are rigorous. 
% For a given interval X and a function X let Y be the computed
% value of f(X). Then f(x) is in Y for all x in X.
% For example
      
x = intval(1e10); format long
sin(x)
          
%%
% Note that the result is rigorous (try sin(2^1000) or similar).
% For timing comparison consider

format short
n=10000; x=random(n,1); X=intval(x);
tic, exp(x); tapprox = toc
tic, exp(X); trigorous = toc
ratio = trigorous/tapprox
   
%% Complex interval standard functions
% Complex interval standard functions are rigorous as well, 
% for example
      
format long
Z = midrad(3+4i,1e-7); 
R = sin(Z)

%%
% It is mathematically correct, that sin(z) is an element of R for
% every complex z with distance less than or equal to 1e-7 from 3+4i.
      
%% Standard functions with argument out of range
% When entering a real argument leading to a complex value of a 
% standard function, there are three possibilities to be specified
% by intvalinit:')

intvalinit('RealStdFctsExcptnNan'); sqrt(intval(-2))
intvalinit('RealStdFctsExcptnWarn'); sqrt(intval(-2))
intvalinit('RealStdFctsExcptnAuto'); sqrt(intval(-2))   

%% Standard functions with argument out of range  and Brouwer's fixed point theorem
% There is a fourth possibility, which is useful in some applications, that is
% to ignore input arguments out of range.
% Note, however, that in this case
% further usage of a result may lead to incorrect conclusions, for example when 
% applying Brouwer's fixed point theorem. 
%
% Consider f(x)=sqrt(x)-1. This function has no real fixed point. However

f = inline('sqrt(x)-1')
X = infsup(-4,2)
intvalinit('RealStdFctsExcptnIgnore'); 
Y = f(X)

%%
% the interval X = [-4,2] is seemingly mapped by f into itself. To avoid such a
% wrong conclusion, one can check whether an input out of range occurred in
% previous computations:

intvalinit('RealStdFctsExcptnOccurred')

%%
% The flag is reset after checking.

%% A common misuse of interval arithmetic
% The dependency problem is the most serious problem of (naive) interval arithmetic.
% The following procedure:
% 
% " Take some numerical algorithm and replace every operation by its corresponding
% interval operation. Then the computed interval result(s) definitely contain the
% true result which would be obtained without the presence of rounding errors. "
%
% will most certainly fail in practice. Although a true statement (if no exception
% like divide by a zero interval occurs), the computed result interval(s) will, for
% very modest problem size, most certainly be of huge diameter and useless. 
%
% Consider, for example, the triangular matrix T where all elements on and below the
% diagonal are equal to 1, and take a randomly generated right hand side. 
% The following lines do this for dimension n=50:

n = 50;
T = tril(ones(n));
b = randn(n,1);

%%
% Then perform a standard forward substitution to compute an inclusion T\b.
% Note that X is defined to be an interval vector, so all operations are
% indeed interval operations (see above section "Invoking interval operations").

X = intval(zeros(n,1));
for i=1:n
  X(i) = b(i) - T(i,1:i-1)*X(1:i-1);
end
X

%%
% The result is displayed with uncertainty perfectly making visible the loss of accuracy. 
% This is due to one of the most common misuses of interval arithmetic, also 
% called "naive interval arithmetic". For more details and examples cf.
%
%  S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic.
%    Acta Numerica, 19:287-449, 2010. 
%
% to be downloaded from "www.ti3.tuhh.de/rump". Note that the linear
% system is very well-conditioned:

cond(T)

%%
% By the well-known rule of thumb of numerical analysis we expect at least
% 14 correct digits in a floating-point approximation T\b. Using a proper
% (non-naive) method, an inclusion of this quality is indeed achieved:

verifylss(T,b)

%%
% Such methods are called "self-validating methods" or "verification methods".
% For some details see the reference above or
%
% S.M. Rump: Self-validating methods. Linear Algebra and its Applications (LAA), 
% 324:3-13, 2001. 
%
% Due to an improved evaluation of the residual (default option "intvalinit('ImprovedResidual')" ,
% see also function "lssresidual.m") 
% 15 correct decimal digits of the result are computed. 


%% Rigorous solution of linear systems
% The INTLAB linear system solver can be called with "\" or "verifylss".' For example, 
% [bare with me, I am often in Japan where the backslash appears like japanese Yen.] 

n = 100; 
A = randn(n); 
b = A*ones(n,1); 
X = verifylss(A,b);

%%
% generates and solves a randomly generated 100x100 linear system. The inclusion
% and its quality is checked by 
     
X([1:3 98:100])
max( X.rad ./ abs(X.mid) )
      
%%
% which calculates the maximum relative error of the inclusion radius with
% respect to the midpoint. The same is done by
     
max(relerr(X))

%% Accuracy of rigorous linear system solving: Hilbert matrices
% For estimating accuracy, try

format long e
n = 10; 
H = hilb(n); 
b = ones(n,1); 
X = verifylss(H,b)
    
%%
% The notoriously ill-conditioned Hilbert matrix is given by H_ij := 1/(i+j-1). 
% To estimate the accuracy, we use the symbolic toolbox to compute the perturbation
% of the solution when perturbing only the (7,7)-element of the input matrix by 2^(-52):

Hs = sym(H,'f'); 
Hs(7,7) = Hs(7,7)*(1+sym(2^(-52))); 
double( Hs \ b )

%%
% The statement "sym(H,'f')" makes sure that no conversion error appears
% when changing H into symbolic format.
% This tiny perturbation already changes the solution in the fourth place;
% thus the computed inclusion is very accurate.

%% Extremely ill-conditioned linear systems
% By default, all computations in INTLAB are, like in Matlab, performed in double precision. This
% restricts treatable linear systems to a maximum condition number of roughly 10^16. 
%
% Starting with INTLAB Version 7, I rewrote my linear system solver completely. 
% Now, although only double precision is used, linear systems with larger condition numbers are 
% solvable. Consider

format long _
n = 20; A = invhilb(n); 
condA = norm(double(inv(sym(A))))*norm(A)

%%
% The common rule of thumb tells that for a condition number 10^k, an algorithm in double precision 
% should produce 16-k correct digits. In our case this means roughly 16-27=-11 "correct" digits, namely
% none. For a random right hand side Matlab computes

b = A*randn(n,1);
x = A\b

%%
% A corresponding warning indicates the difficulty of the problem. Note that in this case the Matlab
% guess of the condition number is pretty good. 
%
% An verified inclusion of the solution is computed by

X = verifylss(A,b,'illco')

%%
% As expected the Matlab approximation differs significantly from the true values, for some 
% components the sign is incorrect. The maximum relative error of the components of the computed 
% inclusion is

max(relerr(X))

%%
% so that each component is accurately included with at least 12 correct figures.
      
%% Structured linear systems
% In general, intervals are treated as independent quantities. If, however, there are dependencies,
% then taking them into account may shrink the solution set significantly. An example is

   format short
   n = 4;  e = 1e-3; intvalinit('displayinfsup');
   A = midrad( toeplitz([0 1+e e 1+e]),1e-4);
   b = A.mid*ones(n,1);
   Amid = A.mid
   X1 = verifylss(A,b)
   
%%   
% First the matrix has been treated as a general interval matrix without dependencies. 
% Recall that only the midpoint is displayed above; all entries of the interval matrix have a uniform
% tolerance of 1e-4.
%
% Several structural information may be applied to the input matrix, for example,

   X2 = structlss(structure(A,'symmetric'),b);
   X3 = structlss(structure(A,'symmetricToeplitz'),b);
   X4 = structlss(structure(A,'generalToeplitz'),b);
   X5 = structlss(structure(A,'persymmetric'),b);
   X6 = structlss(structure(A,'circulant'),b);
   res = [ X1 X2 X3 X4 X5 X6 ];
   rad(res)

%%
% Here only the radii of the inclusions are displayed. Note that the inclusion may become much
% narrower, in particular treating the input data as a circulant matrix.
   

%% Sparse linear systems
% The following generates a random sparse system with about 9 nonzero elements per row.
  
format short      
n=10000; A=sprand(n,n,2/n)+speye(n); A=A*A'; b = ones(n,1);
      
%%
% The linear system is generated to be symmetric positive definite.
% Before calling the verified linear system solver, the fill-in should
% reduced. The original matrix looks like

p = symamd(A); 
spy(A)
title('sparsity pattern of A')

%%
% whereas after minimum degree reordering the matrix looks like

spy(A(p,p))
title('sparsity pattern of renumbered A')
      
%% 
% The timing for the built-in (floating point) solver compared to the 
% verified solver is as follows:

tic, x = A(p,p)\b(p); toc
      
%%

tic, X = verifylss(A(p,p),b(p)); toc
     
%% Inclusion of eigenvalues and eigenvectors
% To compute verified inclusions of eigenvalue/eigenvector pairs of
% simple or multiple eigenvalues, 
% consider, for example, the famous Wilkinson(21) matrix:
     
format long
W = wilkinson(21);              % generation of the matrix
[V,D] = eig(W);                 % eigenvalue/eigenvector approximations
for k=18:21
  [L,X] = verifyeig(W,D(k,k),V(:,k))        % inclusions for the small eigenvalues
end
      
%% Eigenvalue pairs and invariant subspaces
% The smallest eigenvlues are  10.74619418290332 and 10.74619418290339, where all
% displayed digits are verified to be correct.
% Invariant subspaces of nearby eigenvalues are in general ill-conditioned.
% Nearby eigenvalues can also be regarded as clusters. From the inclusions above
% we can judge how narrow the eigenvalues are. So one of the approximations can
% be used as an approximation of the pair. 
   
[L,X] = verifyeig(W,D(18,18),V(:,18:19))    % inclusion of the 18/19 eigenvalue pair
[L,X] = verifyeig(W,D(20,20),V(:,20:21))    % inclusion of the 20/21 eigenvalue pair

%%
% Note that interval output with uncertainty ("_") is used, so all displayed decimal
% places of the bases of the invariant subspaces are verified to be correct. 
% As explained in section "Output formats of intervals III", the inclusion  
% 10.7461941829034_  of the two smallest eigenvlues reads [10.7461941829033,10.7461941829035],
% thus including the true eigenvalues as displayed above. 
% 
% The mathematical statement is that the displayed intervals for the cluster contain
% (at least) two eigenvalues of the Wilkinson matrix W. The size of the cluster
% is determined by the number of columns of the invariant subspace approximation.

%% Eigenvalues of structured matrices
% As for linear systems, the interval input matrix may be structured. Taking into account such
% structure information may shrink the inclusion. As an example consider

   format short
   e = 1e-3;
   A = midrad( toeplitz([0 1+e -e/2 1+e]),1e-4); 
   [v,d] = eig(A.mid); xs = v(:,2:3); lambda = d(2,2);
   X1 = verifyeig(A,lambda,xs);
   X2 = structeig(structure(A,'symmetric'),lambda,xs);
   X3 = structeig(structure(A,'symmetricToeplitz'),lambda,xs);
   X4 = structeig(structure(A,'generalToeplitz'),lambda,xs);
   X5 = structeig(structure(A,'persymmetric'),lambda,xs);
   X6 = structeig(structure(A,'circulant'),lambda,xs);
   res = [ X1 X2 X3 X4 X5 X6 ]; 
   rad(res)
   
%%
% As for linear systems, only the radii of the inclusions are displayed. 


%% Nonlinear systems of equations, polynomials, etc.
% For inclusions of systems of nonlinear equations, of roots of polynomials etc.
% cf. the corresponding demos.

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de
