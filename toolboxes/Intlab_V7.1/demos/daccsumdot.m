%% DEMOACCSUMDOT  A little demonstration of accurate summation and dot products

%% Accurate summation
% Recently I published various algorithms for the accurate computation of sums and dot products.
% A typical application is the accurate computation of a residual.
%
% First we define an ill-conditioned matrix, the inverse Hilbert matrix as provided by Matlab
% (we display only the upper left corner of A):

format short
n = 10;
A = invhilb(n);
v = 1:4;
A(v,v)
cond(A)

%%
% We calculate a right hand side such that the solution is the vector of 1's. Since the matrix
% entries are not too large integers, the true solution is indeed the vector of 1's.
%
% The approximate solution by the built-in Matlab routine is moderately accurate, as expected by
% the condition number of the matrix.

b = A*ones(n,1);
xs = A\b

%% Residual iteration
% If a residual iteration is performed in working precision, the result becomes backward stable;
% however, the forward error does not improve. We display the result after five iterations.

for i=1:5
  xs = xs - A\(A*xs-b);
end
xs

%% Accurate residual iteration
% This changes dramatically when calculating the residual in double the working precision.

format long
xs = xs - A\Dot_(A,xs,-1,b)

%%
% As expected, the accuracy increases. After four iterations the approximation is full accuracy.

for i=1:3
  xs = xs - A\Dot_(A,xs,-1,b);
end
xs

%%
% Note that the residual is calculated "as in" double the working precision, but the result is 
% stored in working precision.

%% Verified inclusion
% The same principle is used in the verification routine verifylss. There is a choice how to
% calculate the residual:

intvalinit('ImprovedResidual')
X1 = verifylss(A,b)

%%
% A heuristic is used to improved the accuracy. It is fast, but not necessarily accurate ("poor
% men's residual"). Calculating the residual as above is slower but more accurate:

intvalinit('QuadrupleResidual')
X2 = verifylss(A,b)

%% Very ill-conditioned matrices
% Next we use an extremely ill-conditioned matrix proposed by Boothroyd (we show some entries
% of the upper left corner). As before the right hand side is computed such that the exact 
% solution is the vector of all 1's.

n = 15;
[A,Ainv] = Boothroyd(n);
A(v,v)
b = A*ones(n,1);

%%
% Since the inverse is the original matrix with a checkerboard sign distribution and thus explicitly
% known, we can safely compute the condition number.

format short
cnd = norm(A)*norm(Ainv)

%%
% As expected, the Matlab approximation has no correct digit, even the sign is not correct.

xs = A\b

%%
% Using accurate dot products based on error-free transformations, an inclusion of the solution can
% be calculated:

format long _
X = verifylss(A,b,'illco')

%% Extremely ill-conditioned sums and dot products
% There are routines to generate extremely ill-conditioned sums and dot products. Consider

n = 50;
cnd = 1e25;
[x,y,c] = GenDot(n,cnd);

%% Computation "as if" in K-fold precision
% Vectors x and y of length n are generated such the condition number of the dot product is cnd=1e25 and the
% true value of x'*y is c. Therefore it can be expected that a floating-point approximation has no
% correct digit, in fact true result and approximation differ significantly in magnitude:

c
x'*y

%%
% The computation of x'*y in double the working precision  gives a more accurate approximation:

c
Dot_(x',y)

%%
% A result "as if" computed in triple the working precision and rounded into working precision is
% accurate to the last bit:

c
Dot_(x',y,3)

%% Accurate approximation
% An alternative is to use error-free transformation to compute an accurate result, independent of
% the condition number. For an extremely ill-conditioned dot product with condition number 1e100 the
% result is still accurate to the last bit.

n = 50;
cnd = 1e100;
[x,y,c] = GenDot(n,cnd);
c
AccDot(x',y)

%%
% An inclusion of the result can be computed as well:

c
infsup(AccDot(x',y,[]))

%% Hidden line
% There is quite some effort in computer geometry to design properly working hidden line algorithms.
% Of course, the decision whether a point is visible or not is simply decided by the sign of some
% dot product. It seems hard to believe, but evaluating dot products in double precision is
% sometimes not enough to make the right decision. In that case an accurate dot product may help.
%
% The following graph shows the solution set of an interval linear system as on the cover of
% Arnold's book. When executing this in Matlab and rotating the graph, sometimes the display is not
% correct.

format short
A = ones(3)*infsup(0,2); A(1:4:end) = 3.5
b = ones(3,1)*infsup(-1,1)
plotlinsol(A,b)
view(-200,20)

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de

 