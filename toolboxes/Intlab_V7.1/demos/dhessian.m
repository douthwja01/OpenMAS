%% DEMOHESSIAN  Short demonstration of Hessians
%

%% Some sample applications of the Hessian toolbox
% Hessians implement second order automatic differentiation in forward mode, which is
% conveniently to implement using the Matlab operator concept.
%

format compact long infsup
setround(0)                           % set rounding to nearest

%% Initialization of Hessians
% The initialization follows the same principles as for gradients, for example

x = hessianinit([ -.3 ; pi ])

%% Access of the Hessian
% Define the function f: R^n->R^n by

f = @(x)( 3*x(1)*x + sin(x).*(sec(x)-sqrt(x)) )

%%
% The number of unknowns is determined by the length of the input vector x.
% For example,

f(1:4)

%%
% The function can be evaluated using the Hessian package as follows:

y = f(hessianinit([1.1 -2.7 pi]'));

%%
% There is direct access of the Hessian of y=f(x) by

y.hx

%%
% However, y.hx contains the Hessians of all three component functions of
% the original function f. To access the Hessians it is recommended to
% use the Hessian of individual components, not components of y.hx:

H3 = y(3).hx

%%
% The matrix H3 is the Hessian of the third component function of f at x.

%% A simple example
% Consider the following function
 
f = inline(' sin(x-log(x+2))-x*cosh(x)/15 ')
 
%% 
% To plot the function first vectorize f :
 
F = vectorize(f)
 
%% 
% Plot the function including the x-axis:
 
x = linspace(-1,3); close; plot( x,F(x), x,0*x )
 
%% Zeros of a function
% There are two roots near 1.5 and 2.5, and two extrema 
% near -0.5 and 2. The roots can be included by verifynlss.
% For this simple function the inclusion is of maximum accuracy.
 
X1 = verifynlss(f,1.5)
X2 = verifynlss(f,2.5)
  
%% Extrema of a function
% The extrema can be approximated and included using Hessians.
% The following is one step of a simple Newton iteration starting at x=-0.5 :
 
x = -0.5; 
y = f(hessianinit(x)); 
x = x - y.hx\y.dx';
y
  
%% Inclusion of extrema
% Inclusions of the extrema of f are computed by "verifynlss" with 
% parameter 'h': This parameter specifies that
% f'(x) = 0 is solved instead of f(x) = 0.
 
X1 = verifynlss(f,-0.5,'h')
X2 = verifynlss(f,2,'h')
 
%% Functions in several unknowns
% Function with several unknowns are handled in a similar manner. 
% Consider the following test function by N. Gould. It is taken from
% the Coconut collection of test function for global optimization, 
% http://www.mat.univie.ac.at/~neum/glopt/coconut/benchmark/Library2.html .
 
f = inline(' x(3)-1 + sqr(x(1)) + sqr(x(2)) + sqr(x(3)+x(4)) + sqr(sin(x(3))) +  sqr(x(1))*sqr(x(2)) + x(4)-3 + sqr(sin(x(3))) + sqr(x(4)-1) + sqr(sqr(x(2))) + sqr(sqr(x(3)) + sqr(x(4)+x(1))) + sqr(x(1)-4 + sqr(sin(x(4))) + sqr(x(2))*sqr(x(3))) + sqr(sqr(sin(x(4)))) ')
 
%% Approximation of an extremum
% On the Web-site the global minimum of that function in 4 unknowns is 
% given as 5.7444 . We use a couple of a Newton iterations 
% starting at x=ones(4,1) to approximate a stationary point: 

format short
x = ones(4,1); 
for i=1:18
   y = f(hessianinit(x)); 
   x = x - y.hx\y.dx';
end
y.dx

%%
% Now x is an approximation of a stationary point:
% The gradient of f evaluated at x is not too far from zero.
   
%% Inclusion of a stationary point
% Using this approximation an inclusion of a stationary point of f is 
% computed by (in this case the last parameter 1 is used to see 
% intermediate results):
 
format long
X = verifynlss(f,x,'h',1)
 
%% Inclusion of a minimum
% The interval vector X includes a root of f', i.e. a stationary point xx of f. To prove
% that f has a minumum at xx we need to prove positive definiteness of
% the Hessian of f evaluated at xx. The interval vector X includes the stationary point
% xx of f. So we compute an inclusion Y of the Hessian at X. 
%
% Mathematically,
% for every x in X the following is true: Y.x is an inclusion of f(x),
% Y.dx is an inclusion of f'(x), and Y.hx is an inclusion of the
% Hessian of f at x. Especially, the Hessian of f evaluated at xx is
% included in Y.hx. 
 
Y = f(hessianinit(X)); 
format _
Y.hx
  
%% 
% This interval matrix contains obviously only diagonally dominant
% matrices, so the stationary point xx of f in X is indeed a (local)
% minimum.  
%

%% A model problem in 5000 unknowns I
% The next problem is taken from 
%
% http://www.sor.princeton.edu/~rvdb/ampl/nlmodels/cute/bdqrtic.mod
%
%    Source: Problem 61 in
%       A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
%       "Performance of a multifrontal scheme for partially separable optimization",
%        Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
%        Copyright (C) 2001 Princeton University
%        All Rights Reserved
%    see bottom of file test_h.m
% 
% The model problem is
%  
%     N = length(x);      % model problem: N = 1000, initial approximation x=ones(N,1);
%     I = 1:N-4;
%     y = sum( (-4*x(I)+3.0).^2 ) + sum( ( x(I).^2 + 2*x(I+1).^2 + ...
%               3*x(I+2).^2 + 4*x(I+3).^2 + 5*x(N).^2 ).^2 );
% 
% This function is evaluated by
%
%     index = 2;
%     y = test_h(x,index);

%% A model problem in 5000 unknowns II
% The given starting vector is x = ones(5000,1).
% Recall that y = f(hessianinit(x)) has 5000 elements in y.x, 2.5e7 elements
% in y.dx and 1.25e11 elements in y.hx. In full storage this would mean 1 TeraByte
% of storage. 
%
% The problem can be solved in the above manner. On my 2.8 GHz
% Laptop this requires about 1 second per Newton iteration, and 5 seconds
% for one function evaluation with interval argument. 
%
% The following calculates an inclusion of a stationary point of f by first
% performing a simple Newton iteration followed by a verification step for
% the nonlinear system. The Hessian is treated as full matrix, so the 
% computation may take a while.
 
n = 5000;
index = 2;
tic
X = verifynlss('test_h',ones(n,1),'h',1,index);
tfull = toc
max(relerr(X))
         
%% A model problem in 5000 unknowns III
% Note the small maximum relative error of the inclusion of the result.
% Verification is faster when solving the nonlinear system treating the
% Hessian as sparse. This is done by the following. 
 
n = 5000;
index = 2;
tic
X = verifynlss('test_h',ones(n,1),'hSparseSPD',1,index);
tsparse = toc
tfull
max(relerr(X))
       
%% 
% Note that verification is faster at the price of a less narrow inclusion
% (for comparison the previous time tfull is displayed). 
% 

%% Verification of a minimum
% Finally, the Hessian at X is computed which includes the Hessian at the
% stationary point xhat in X. 
 
y = test_h(hessianinit(X),index); 
isspd(y.hx)
      
%%
% The interval Hessian is proved by "isspd" to be symmetric
% positive definite: Mathematically the result "isspd(M)=1" for a symmetric
% (Hermitian) interval matrix M proves that every symmetric (Hermitian) matrix
% A within M is positive definite. In our case especially the Hermitian of f
% at the stationary point xhat. Therefore, xhat is proved to be a local minimum of f.

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de
