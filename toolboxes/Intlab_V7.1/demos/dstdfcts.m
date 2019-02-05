%% DEMOSTDFCTS  Interval standard functions in INTLAB
%

%%
% As explained in "DEMOINTVAL" standard functions with interval argument
% compute an inclusion of the true function value or the true function
% range.

%% Input arguments I
% 
% Special care is necessary if an input argument is not a floating-point
% number. For example,

format long _
x = intval([2.5 -0.3]);
y = erf(x)

%%
% computes the true range of the error function erf(2.5), but not
% necessarily of erf(-0.3). To obtain true results use

x = intval('2.5 -0.3');
y = erf(x)

%%
% Note that intval conversts a character string representing a vector
% always into a column vector. 
% To see the accuracy of the result, the infsup-notation may be used:

infsup(y)

%% Input arguments II
% Similarly, the range of a function is included by

x = intval('[2.5,2.500001]')
y = erf(x)

%%
% Note that the large diameter of the output is due to the large diameter
% of the input.
      
%% Accuracy of standard functions
% The complementary error function erfc(x) is defined by 1-erf(x). It
% rapidly approximates 1 for larger x. In this case erfc(x) is much more
% accurate than 1-erf(x).

x = intval(5);
y1 = erfc(x)
y2 = 1 - erf(x)
infsup([y1;y2])

%% Nonlinear equations involving standard functions I
% The zero of a nonlinear function may be approximated by some Newton
% procedure. Consider

x = linspace(0,10);
f = inline('erf(x)-sin(x)')
close
plot(x,f(x),x,0);

%%
% The graph shows the function erf(x)-sin(x) between 0 and 10. Besides the
% obvious root zero there are two small roots near 1 and 2, and
% mathematically it follows that there must be two roots near 8. 
% To approximate the roots, a Newton procedure may be used. For example

x = 1;
for i=1:5
  y = f(gradientinit(x));
  x = x - y.x/y.dx
end

%% Inclusion of the solution of nonlinear systems with standard functions
% As expected the iteration converges very rapidly. The final value is a
% very good approximation of a root of f, and an inclusion is computed by

Y = verifynlss(f,x)

%% Nonlinear equations involving standard functions II
% More interesting is the root cluster near x=8. An attempt to compute
% inclusions fails:

Y1 = verifynlss(f,7.5)
Y2 = verifynlss(f,8.5)

%%
% The reason is that numerically this is a double zero because erf(8) is
% very close to 1, as can be checked by erfc(8) = 1-erf(8) :

y = erfc(intval(8))

%%
% At least an inclusion of the extremum of f near 8 can be computed.
% Clearly it must be near to 2.5*pi.

Y = verifynlss(f,8,'h')
y = 2.5*pi

%%
% For interval enthusiasts an inclusion of 2.5*pi can be computed as well.
% Note, however, that this is only close to the extremum.

y1 = 2.5*intval('pi')

%% Accuracy of the gamma function I
% Usually Matlab standard functions produce very accurate approximations.
% For example,

x = pi
y = gamma(pi)
Y = gamma(intval('pi'))

%%
% The latter Y is a true inclusion of the value of the gamma function at
% the transcendental number pi: X=intval('pi') is an inclusion of the true
% pi, and the gamma function for an interval argument produces an inclusion
% of all values, in particular for the transcendental number pi.

%% Accuracy of the gamma function II
% However, for negative values near integers, where the gamma function has
% poles, Matlab's built-in gamma function is sometimes very inaccurate. For
% some reason this happens only right to negative integers. For example,

e = 2^(-40);
x = [-1-e -1+e];
gamma(x)
gamma(intval(x))

%%
% In this case the input arguments -1+/-e are exactly representable
% floating-point numers. For the example this is not important because both
% the Matlab's approximate gamma function and INTLAB's verified gamma
% function use the same input.
%
% As can be seen only 4 digits are correct for the second argument. Getting
% closer to the right of negative integers makes things worse:

x = succ(-1);
gamma(x)
gamma(intval(x))

%%
% Now the approximate value has no correct digit. If a symbolic toolbox is
% available, the accuracy of INTLAB's gamma function can be checked:

digits(100)
vpa( gamma(sym(x,'f')) )

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de
