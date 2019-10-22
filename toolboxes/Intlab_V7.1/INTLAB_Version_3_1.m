%Help file for INTLAB Version 3.1
%
%
%Additions/changes:
%
%  - slopes for matrices supported
%  - rounding unchanged after use of INTLAB routines
%  - midradcmplx replaced by cintval
%  - correct setround .m/.dll routine detected automatically
%  - results of comparison ==, <, etc. of same size as operands (like in Matlab)
%
%  Others:
%
%  - various performance improvements
%  - certain utility routines added
%  - displaywidth to change width of display
%  - function disp2str added
%  - message 'no inclusion achieved' omitted (output NaN anyway)
%  - For safety, INTLAB terminates if rounding does not work properly
%  - Display window for generating std fcts omitted
%
%
%
%
%We added a nonlinear demo function, showing ease of use but also
% interpretation overhead. In
%   J. More: A Collection of Nonlinear Model Problems. In: Computational
%     Solution of Nonlinear Systems of Equations (Eds.: E.L. Allgower,
%     K. Georg), Lectures in Applied Mathematics, Volume 26, American
%     Mathematical Society, Providence, Rhode Island (1990)
%Fletcher describes a "distillation column test problem". Fletcher gives
% data for three processes, for the third (the methanol-8 problem with
% 31 unknowns) he writes: "I still do not know if there exists a solution
% to this problem."
% The best he found was an approximation with 1e-2 residual norm.
%
%All three problem where solved with verified bounds in
%   G. Alefeld, A. Gienger and F. Potra: Efficient Numerical Validation
%     of Solutions of Nonlinear Systems, SIAM J. Num. Anal. 31, 252-260 (1994).
%
%The definition of Fletcher's test problems is a little involved.
% We added a corresponding program to show the ease of use of INTLAB.
% The call
%    xs = fletcher('init1')
%generates the initial approximation for the first problem (hydrocarbon-6,
% 29 unknowns), whereas input arguments 'init2' and 'init3' generate the
% initial approximation for the second (hydrocarbon-20, 99 unknowns) and the
% third (methanol-8, 31 unknowns) problem, respectively.
%All three problem where solved with verified bounds. The call is simply
%    X = verifynlss('fletcher',xs)
%The number of the problem is choosen in fletcher.m according the dimension
% of xs.
%
%                                  rel. error
%    problem         unknowns   median  maximum   time (300 MHz CPU)
%=======================================================================
% 1  hydrocarbon-6      29     9.2e-14  3.6e-13        2.5 sec
% 2  hydrocarbon-20     99     2.8e-12  8.4e-11         19 sec
% 3  methanol-8         31     3.8e-14  2.7e-13        3.0 sec
%
%The computing time is measured with my a little aged 750 MHz Laptop.
% Switching from gradient to slope evaluation by
%    X = verifynlss('fletcher',xs,'s')
%does not improve the result but takes three times as long. Interpretation
% overhead for nonlinear problems can be quite severe.
%
%
%Finally another note to ease of use. I was asked, whether it is not possible to
%solve nonlinear systems such as Broyden's example
%
% .5*sin(x1*x2) - x2/(4*pi) - x1/2  =  0
%  (1-1/(4*pi))*(exp(2*x1)-exp(1)) + exp(1)*x2/pi - 2*exp(1)*x1 )  =  0
%
%with initial approximation [ .6 ; 3 ] and one solution [ .5 ; pi ] without use of
%a file. Yes, it is. Define
%
%  str = '[ .5*sin(x(1)*x(2)) - x(2)/(4*P1) - x(1)/2 ; ...
%              (1-1/(4*P1))*(exp(2*x(1))-exp(1)) + exp(1)*x(2)/P1 - 2*exp(1)*x(1) ]' ;
%  f = inline( str , 1);
%
%This specifies f to be a function in the variable x (which is a vector of dimension 2)
%and one parameter P1, the number pi. Of course, f can be defined in one statement.
%Then call
%
%   format long
%   xs = [ .6 ; 3 ] ;  P1 = midrad(pi,1e-15);
%   Y = verifynlss(f,xs,'g',0,P1)
%
%to obtain
%
%   intval Y =
%      0.50000000000___
%      3.1415926536____
%
%For further details see "inline" and "verifynlss".
