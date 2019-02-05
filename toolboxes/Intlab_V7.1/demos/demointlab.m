function demointlab
%DEMOINTLAB   Wrapper routine to call INTLAB demos
%
%A selection of INTLAB demos, call
%
%  demointlab
%

% written  10/13/12     S.M. Rump  
% modified 11/07/12     S.M. Rump  INTLAB_larger added
%

d = which('demointlab');
addpath([ d(1:end-13) '\html' ])

clc
disp('Welcome to INTLAB, the Matlab toolbox for reliable computing.')
disp(' ')
disp('The current version 7 consists of more than 800 .m-functions with more ')
disp('  than 30 thousand lines of Matlab-code (more than 50 KLOC with comments). ')
disp('The test suite for INTLAB consists of another 60 KLOC. ')
disp(' ')
disp(' ')
while 1
  displaycomments
  str = input('select demo ','s');
  switch lower(str)
    case '1', web('dintlab_larger.html');
    case '2', web('dintlab.html');
    case '3', web('dintval.html');
    case '4', web('darithmetic.html');
    case '5', web('daccsumdot.html');
    case '6', web('dgradient.html');
    case '7', web('dhessian.html');
    case '8', web('dslope.html');
    case '9', web('dtaylor.html');
    case 'a', web('dpolynom.html');
    case 'b', web('dlong.html');
    case '0', break;
  end
end

disp(' ')
disp('Enjoy INTLAB. Comments and suggestions always welcome to rump (at) tuhh.de .')
disp(' ')

  
  
function displaycomments
  disp(' ')
  disp('This is a wrapper routine to call several INTLAB demos, selected by numbers. ')
  disp(' ')
  disp('1  Some larger examples with INTLAB')
  disp('2  A general demo of some features of INTLAB')
  disp('3  Some examples of interval computations')
  disp('4  Details about interval arithmetic')
  disp('5  Accurate summation and dot products')
  disp('6  The gradient toolbox (gradients of multivariate functions)')
  disp('7  The Hessian toolbox (Hessians of multivariate functions)')
  disp('8  The slope toolbox (slope of multivariate functions')
  disp('9  The Taylor toolbox (taylor expansion of univariate functions)')
  disp('a  The polynomial toolbox (univariate and multivariate polynomials')
  disp('b  The long number toolbox (a rudemantary implementation, originally for internal use)')
  disp(' ')
  disp('0  exit this wrapper')
  disp(' ')
