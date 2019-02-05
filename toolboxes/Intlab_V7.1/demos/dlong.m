%% DEMOLONG  Short demonstration of long numbers
%

%%
% The purpose of the long toolbox was to compute rigorous bounds for the
% value of certain standard functions. Those values are needed to initialize
% the INTLAB system. The long toolbox is slow, but fast enough to do the job.

setround(0)                 % set rounding to nearest
longprecision(0);           % default option

%% Definition of long numbers         
%
% Long numbers are stored in midpoint radius representation. 
% The midpoint is stored in an array of precision to be specified,   
% the error is stored in one number. Long numbers or vectors are generated 
% by the constructor long: 

x = long(7)
V = long([-1;3])

%%
% Long vectors are always column vectors (or forced to column vectors). 

%% Conversion 
% Conversion of double to long by the constructor "long" is always into
% the current internal precision.
%  
% The internal precision may be specified by "longprecision". 
% The call without input parameter gives the current working precision 
% in decimals, the call with input parameter sets working precision. 

p = longprecision
longprecision(50)

%% 
% Now, the statement "x = long(7)" generates a long number with  
% approximately 50 decimal digits. This is only approximate because 
% internal representation is to some base beta, a power of 2.

%% Output of long numbers
% Output is usually to a little more than double precision. If you want 
% to see more digits, say k, use "display" with second parameter equal 
% to k. To see all digits, use k=0.

longprecision
x = 1/long(7)
display(x,40)
display(x,0)

%%
% Output of long numbers is not rigorous. All but a few of the last
% digits are correct. 

%% Arithmetic operations
% Long operations +,-,*,/ and ^ are supported. Note that operations on 
% vectors are always performed elementwise.

x = [ long(3) ; -7 ]
x*x

%% Output of long intervals
% The display routine takes uncertainties into account. Only the 
% correct digits plus some extra are displayed.

longprecision(50); 
x = long(1)/37; 
display(x,0)
for i=1:100
  x=x*x; x=x*37; 
end
display(x,0)

%% Interval and non-interval operations
% Computing with uncertainties may be switched off by 

longinit('WithoutErrorTerm')
longprecision(50); 
x = long(1)/37; 
display(x,0)
for i=1:100
  x=x*x; x=x*37; 
end
display(x,0)

%%
% In this case all digits including incorrect ones are displayed.
% Computing without error term is a usual long precision arithmetic
% with specified precision. Note that scalar operations suffer from
% quite some interpretation overhead.

%% Conversion between long and double 
% Conversion from long to double is approximately to nearest, conversion 
% to interval is rigorous. 
%
% For example, in the following the function "longpi" calculates
% "pi" to the specified longprecision, "IntPi" is a true inclusion of the 
% transcendental number "pi".

longinit('WithErrorTerm'); 
longprecision(100); 
Pi = longpi;
display(Pi,0)
flptPi = long2dble(Pi)
IntPi = long2intval(Pi)
format long
infsup(IntPi)

%% Long numbers with error term
% Long numbers may be specified with an explicit error term.
% For example, 

longprecision(50); 
x = long(-1.5)
display(x,0)
x = addlongerror(x,1e-40)
display(x,0)

%%
% defines x to be an interval with midpoint -1.5 and radius 
% approximately 10^(-40). Only meaningful digits are stored and displayed.

%% Specifying extremely small errors
% For very small errors leaving the range double precision 
% floating point numbers, the error may be specified by 
% the mantissa and the exponent of the error:

longprecision(50); 
x = long(2^-1000)^2; 
x = addlongerror(x,1,-620)

%%
% The final x, which is 2^(-2000), is afflicted with an error
% of 10^(-620).

%% Taylor series: an example
% As an example, the following code computes the value of E = exp(x)
% by a Taylor series:

p = 100; longprecision(p); 
x = -30;
t = 1; T = long(1); E = T; k = 0;
while abs(t)>10^(-p)
  k = k+1;
  t = t*x/k;
  T = T*x/k;
  E = E + T;
end
k
exp(x)
display(E,0)

%%
% Note that for large negative values of x there quite some 
% cancellation. This can be seen by

x = 30;
t = 1; T = long(1); E = T; k = 0;
while abs(t)>10^(-p)
  k = k+1;
  t = t*x/k;
  T = T*x/k;
  E = E + T;
end
k
1/exp(x)
display(1/E,0)

%% Ill-conditioned polynomials
% Consider the following polynomial:

P = inline(' 4999*x.^6 - 200*x.^5 + 102*x.^4 - 2*x.^3 - 2500*x.^2 + 100*x - 1 ')

%%
% This is an example Bugeaud-Mignotte polynomial. The general form is
%
% ( X^n - aX + 1 )^k - 2X^(nk-k)(aX-1)^k
%
% where a>=10, n>=3 and k>=2.
%
% Those polynomials are constructed to have a pair of very close real roots near c=1/a+1/a^(n+1). 
% A graph near c looks as follows:

e = 3e-8; 
c = 1/50+1/50^4;
x = c*(1+linspace(-e,e));
close
plot(x,P(x),x,0*x)

%%
% From the graph it is not clear whether the polynomial has no, a double or two real roots
% in the interval c*[1-e,1+e]. An evaluation using the long package yields
% the following:

y = long2dble(P(long(x)));
close
plot(x,y,x,0*x)

%% Sample programs
% For sample programs using long numbers, see for example the
% source codes of long\longpi.m or long\@long\exp.m  .

%% Enjoy INTLAB
% INTLAB was designed and written by S.M. Rump, head of the Institute for Reliable Computing,
% Hamburg University of Technology. Suggestions are always welcome to rump (at) tuhh.de
