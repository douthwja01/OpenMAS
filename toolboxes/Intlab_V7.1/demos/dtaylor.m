%% DEMOTAYLOR    Short demonstration of Taylor toolbox
%

%% Some sample applications of the Taylor toolbox
% The Taylor toolbox computes Taylor coefficients of a univariate function 
% in forward mode, which is conveniently to use by the Matlab operator concept.
% It works much in the spirit of the gradient and Hessian toolbox, so
% I recommend to visit the gradient and Hessian demo first.
%

%% Initialization of Taylor variables
% In order to use the automatic Taylor toolbox, the independent variable need
% to be identified and a value has to be assigned. This is performed by
% the function "taylorinit", for example                            

format compact short _
u = taylorinit(2.89)

%% Operations between Taylor variables
% If at least one operand is of type Taylor, operations are executed as 
% Taylor operations.
% For example,                                                          

x = taylorinit(3.5);  
y = sin(3*x-sqrt(x+5))

%% Interval Taylor variables
% If arguments are of type intval, an inclusion of the true value is computed:

x = taylorinit(midrad(3.5,1e-12));  
y = sin(3*x-sqrt(x+5))

%%
% For f(x):=exp(3*x-sqrt(x)), the result y contains in y.t the function value f(3.5)
% and the first 4 derivatives:
% 

y.t

%% Taylor vector variables
% Note that the Taylor toolbox accepts one independent variable. One may initialize
% a Taylor variable of a vector argument; this is the same as initializing each
% component as the independent variable (with a different value). It is convenient
% for function evaluations with many arguments:

f = inline('sin(3*x-sqrt(x+5))')
x = taylorinit([-3 0.1 3.5]')
y = f(x)


%% Complex arguments
% When evaluating the expression for another argument, use the same
% statement as before with new values. Here we assign the Taylor variable to 
% carry 2 derivatives (the default is 4):

x = taylorinit(-3.5+.2i,4);  
y = sin(3*x-sqrt(x))

%% Access to the Taylor coefficients
% The Taylor coefficients are accessed by {}, so that y{0} is the function value and
% y{k} denotes the k-th Taylor coefficient:

y{0}
y{1:3}

%% Access to the Taylor coefficients of vector variables
% When initializing a Taylor vector, the individual vector components are accessed
% by () and Taylor coefficients by {}. For example,

f = inline('sin(3*x-sqrt(x+5))')
x = taylorinit([-3 0.1 3.5]')
y = f(x)
y(1)
y{2}(3)

%%
% accesses the Taylor value f(-3) and the second Taylor coefficient of f(3.5), 
% respectively. 

%% An example: Taylor series          
% Define 

f = inline('sinh(x-exp(2/x))')

%%
% Then the Taylor coefficients 0..7 of f at x=1.234 are computed by

kmax = 4;
x = 1.234;
y = f(taylorinit(x,kmax))

%%
% The Taylor coefficients y{k} satisfy
%   f(x+e) = sum[0..k]( y{k}*e^k ) + O(e^(k+1)) :

format long
e = 1e-3;
v = f(x+e)
yapprox = sum( y{0:kmax} .* e.^(0:kmax) )

%% Inclusion of function value by Taylor series
% For an inclusion of the function value we may calculate the Taylor coefficients
% in interval arithmetic and add the error term:

format long _
x = intval('1.234');
Y = f(taylorinit(x,kmax));
e = intval('1e-3');
Y_ = f(taylorinit(x+hull(0,e),kmax+1));
for k=0:kmax
  Yincl = sum( Y{0:k} .* e.^(0:k) ) + Y_{k+1}*e.^(k+1)
end

%%
% Note how nicely the linear convergence can be observed by the "_"-notation. Also
% note that this is a true inclusion of f(1.234+1e-3)=f(1.235) because both arguments
% x=1.234 and e=1e-3 are intervals including the decimal numbers 1.234 and 0.001
% (both are not floating-point numbers).

%% An Application: integration
% Consider

f = inline('sin(pi*x)-sin(x)'); a = 0; b = 20;
x = linspace(a,b,1000); close, plot(x,f(x),x,0)

%%
% It is easy to see that for the transcendental number pi, the true value of the 
% integral of f from a to b is cos(b)-1:

cos(b)-1

%%
% There is a rudemtary integration routine "verifyquad" using Romberg's roule 
% based on the Taylor toolbox. It calculates 

ApproxIncl = verifyquad(f,a,b)
infsup(ApproxIncl)

%% Integration using transcendental constants
% This is a true inclusion of the integral with "pi" denoting the floating-point 
% approximation of the transcendental number pi. To calculate an inclusion of
% the function with the true transcendental number pi, we use the following
% program:

% function y = testfuntaylor(x)
%   if isintval(x)
%     Pi = 4*atan(intval(1));
%   else
%     Pi = pi;
%   end
%   y = sin(Pi*x)-sin(x);
%

%%
% The result, however, does not change very much:

TrueIncl = verifyquad(@testfuntaylor,a,b)
infsup(TrueIncl)

%% A comparison to the Matlab function "quad"
% For this particular function the approximate routine may get problems if we specify
% a little more accuracy:

e = 1e-12;
tic, Approx = quad(@testfuntaylor,a,b,e), toc,
tic, Incl = verifyquad(@testfuntaylor,a,b,e), toc

%%
% Note that the verification routine is faster and calculates an inclusion of the
% 'true' function (with the transcendental number pi). Insisting on
% even more accuracy make things worse:

e = 1e-14;
tic, Approx = quad(@testfuntaylor,a,b,e), toc,
tic, Incl = verifyquad(@testfuntaylor,a,b,e), toc

%%
% Note that the Matlab routine "quad" gives no error message.

%%
% Now the approximate value has no correct digit and the verification routine is
% still faster. However, it may be the other
% way around, at least concerning computing time:

f = inline('sqrt(x)'); a = 0.0001; b = 2;
tic, Approx = quad(f,a,b), toc,
tic, Incl = verifyquad(f,a,b), toc

%%
% Since "verifyquad" requires differentiability properties of the integrand,
% it takes a lot of computing time near a singularity. 

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de
