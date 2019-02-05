function [c,err] = intval(x,y,type)
%INTVAL       Intval class constructor
%
%  c = intval(x)        type cast for x double or intval
%  c = intval(s)        verified conversion for string s
%
%For verified conversion, string s may be one- or two-dimensional
%  containing one or more input numbers. Input may be real or complex,
%  scalars or intervals in
%     infsup representation  [ ... , ... ]  or
%     midrad representation  < ... , ... >  or
%     numbers with uncertainty  3.14_+2.71_i
%Output c is always interval column vector
%
%  Examples:  c = intval('0.1');    ==>  rad(c) = 1.3878e-017
%
%             c = intval('0.1e1');  ==>  rad(c) = 0
%
%             c = intval('3.14_');  ==>  rad(c) = 1e-2
%
%             input   >> s = [ '-3.4_e2 0 ' ; '.02 +123  ' ], intval(s)
%             output  s =
%                     -3.4_e2 0
%                     .02 +123
%                     intval ans =
%                      -34_.____
%                         0.0000
%                         0.0200
%                       123.0000
%
%             input   >> intval(' [1.993,2] <3+2i,1e-3> -4.33_i ')
%             output  intval ans =
%                         2.00__ +  0.00__i
%                         3.0000 +  2.00__i
%                         0.0___ -  4.3___i
%
%The last interval has thick real part because complex intervals are always
%  midpoint/radius intervals, i.e. the uncertainty in the imaginary part
%  -4.33_i is interpreted as a radius in the complex plane.
%
%Certain mathematical constants are predefined. Use
%
%    C = intval(char)       narrowest interval C containing constant
%or
%    [c,C] = intval(char)   fl-pt c and narrowest interval C s.t. const in c+C
%
%for
%
%char   'e'         base of natural logarithm   2.71828182845904523536 ...
%       'pi'        pi                          3.14159265358979323846 ...
%       'c'         Euler-Mascheroni constant   0.57721566490153286060 ...
%       'Conway'    Conway constant             1.30357726903429639125 ...
%       'Catalan'   Catalan constant            0.91596559417721901505 ...
%       'golden'    golden ratio                1.61803398874989484820 ...
%       'log2'      natural logarithm of 2      0.92457068871300903284 ...
%       'log10'     natural logarithm of 10     2.30258509299404568401 ...
%       'sqrt2'     square root of 2            1.41421356237309504880 ...
%       'sqrt3'     square root of 3            1.73205080756887729352 ...
%       'sqrt5'     square root of 5            2.23606797749978969640 ...
%
%
%To avoid error messages use, for example, 
%
%  [x,err] = intval('[3,2]')
%
%where  err =  0   normal end, no error detected
%              1   improper interval (lower bound greater than upper)
%              2   other errors
%
%If one of the operands of a unary or binary arithmetic operation or a
%  standard function is of type intval, the operation is executed with
%  interval arithmetic producing rigorous bounds for the result.
%The Matlab system parses - as every compiler - expressions following the
%  priority rules of the operation. For example,
%
%  x = intval(0.125) + 1/10
%
%does not contain the decimal number 0.225 because 1/10 is executed in
%  floating-point with round to nearest (in fact, x is a point interval,
%  try x.sup-x.inf). A correct answer is produced by
%
%  x = 0.125 + intval(1)/10 .
%
%The easiest way always to produce correct answers is that all quantities
%  occuring in an expression are of type intval.
%
%The result of standard functions with interval arguments is also correct,
%  even for extreme examples like
%
%  sin(intval(2^1000)) .
%
%For more information try demointval. For information of system variables,
%  default output format, settings for standard functions and others try
%
%  help intvalinit 
%
%and the various demos.
%


%
%Internal representation of real intervals by infimum and supremum:
%  [ inf , sup ]  represents the set  { x in R | inf <= x <= sup }
%
%Internal representation of complex intervals by midpoint and radius:
%  < mid , rad >  represents the set  { x in C | abs(x-mid) <= rad }
%
%  inf/sup   may be a real number, vector, matrix, or sparse matrix
%  mid       may be a complex number, vector, matrix, or sparse matrix
%  rad       may be a nonnegative real number, vector, matrix, or sparse matrix
%
%Note that for an array x=1:6, for example, the assignment x(3)=[] results in a
%  vector [1 2 4 5 6] with 5 components. An interval with value NaN or component
%  NaN does not indicate an empty interval but rather that no information about
%  this interval or component is known.
%

%If .complex = 0   real interval by inf/sup, mid=rad=[]
%   .complex = 1   complex interval by mid/rad, inf=sup=[]
%

%Function intval with more than one argument only for internal use.
%To produce a non-point interval use ***only*** functions infsup, midrad or
%  intval with input string
%

% written  10/16/98     S.M. Rump
% modified 12/30/98     S.M. Rump  slopes and bug fix for save/load
% modified 10/24/99     S.M. Rump  sparse matrices, comment empty intervals
% modified 08/29/00     S.M. Rump  multivariate polynom handling
% modified 12/18/02     S.M. Rump  hessians added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    error message for invalid application of 'infsup'
%                                    Matlab sparse bug
% modified 01/01/05     S.M. Rump  error parameter added to intval('..') (thanks to Arnold)
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 09/10/07     S.M. Rump  performance, huge arrays
% modified 10/20/07     S.M. Rump  sparse Matlab 6.5 bug
% modified 06/23/07     S.M. Rump  logical
% modified 11/20/09     S.M. Rump  interval constants
% modified 03/07/10     S.M. Rump  rndold
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==1
    if isa(x,'logical')
      x = double(x);
    end
    if isa(x,'double')          % 1 parameter, point interval
      if isreal(x)
        c.complex = 0;
        c.inf = x;
        c.sup = x;
        [m,n] = size(x);        % take care of thin inf-interval
        c.mid = [];
        c.rad = [];
      else
        c.complex = 1;
        c.inf = [];
        c.sup = [];
        c.mid = x;
        % careful: just 0 is not sparse and may cause tremendous memory need
        if issparse(x)
          c.rad = sparse(size(x,1),size(x,2),0);
        else
          c.rad = 0;
        end
        [m,n] = size(x);        % take care of thin inf-interval
        if m*n<2^31             % not huge input
          index = isinf(x);
          if any(index(:))
            c.mid(index) = 0;
            if isequal(c.rad,0) % c.rad is already array if c.rad is sparse
              c.rad = zeros(size(x));
            end
            c.rad(index) = inf;
          end
        else
          [I,J,S] = find(x);
          index = isinf(S);
          if any(index)         % S is vector, so is index
            S_ = S;
            S_(index) = 0;
            c.mid = sparse(I,J,S_,m,n);
            S(~index) = 0;
            S(index) = inf;     % take care of -inf
            c.rad = sparse(I,J,S,m,n);
          end
        end
      end

      c = class(c,'intval');
      
    elseif isa(x,'intval')      % 1 parameter, interval
      c = x;

    elseif isa(x,'char')        % 1 parameter input string, verified conversion
      
      switch lower(x)           % predefined math constants
        case 'e'
          % 2.718281828459045235360287471352662497757247093699959574966967 ...
          if nargout<=1
            Cinf = 91210402*2^(-25)+103329490*2^(-52);
            Csup = 91210402*2^(-25)+103329492*2^(-52);
          else
            C = 91210402*2^(-25)+103329490*2^(-52);
            Cinf = 87383992*2^(-79)+90309234*2^(-106);
            Csup = 87383992*2^(-79)+90309236*2^(-106);
          end
        case 'pi'
          % 3.1415926535897932384626433832795028841971693993751058209749445 ...
          if nargout<=1
            Cinf = 105414357*2^(-25)+71487872*2^(-55);
            Csup = 105414357*2^(-25)+71487888*2^(-55);
          else
            C = 105414357*2^(-25)+71487872*2^(-55);
            Cinf = 74025356*2^(-79)+103331852*2^(-106);
            Csup = 74025356*2^(-79)+103331854*2^(-106);
          end
        case 'c'
          % 0.57721566490153286060651209008240243104215933593992 ...
          if nargout<=1
            Cinf = 77472575*2^(-27)+117137792*2^(-57);
            Csup = 77472575*2^(-27)+117137808*2^(-57);
          else
            C = 77472575*2^(-27)+117137808*2^(-57);
            Cinf = -95609884*2^(-84)-133163776*2^(-116);
            Csup = -95609884*2^(-84)-133163712*2^(-116);
          end
        case 'conway'
          % 1.3035772690342963912570991121525518907307025046594048757548613906285
          if nargout<=1
            Cinf = 87481589*2^(-26)+88733220*2^(-53);
            Csup = 87481589*2^(-26)+88733222*2^(-53);
          else
            C = 87481589*2^(-26)+88733220*2^(-53);
            Cinf = 89638218*2^(-83)+90076514*2^(-110);
            Csup = 89638218*2^(-83)+90076516*2^(-110);
          end
        case 'catalan'
          % 0.915965594177219015054603514932384110774 ...
          if nargout<=1
            Cinf = 122938820*2^(-27)+131081914*2^(-54);
            Csup = 122938820*2^(-27)+131081916*2^(-54);
          else
            C = 122938820*2^(-27)+131081914*2^(-54);
            Cinf = 72488322*2^(-84)+96926304*2^(-113);
            Csup = 72488322*2^(-84)+96926312*2^(-113);
          end
        case 'golden'
          % 1.6180339887498948482045868343656381177203091798057628621354486227052604 ...
          if nargout<=1
            Cinf = 108584422*2^(-26)+120580430*2^(-53);
            Csup = 108584422*2^(-26)+120580432*2^(-53);
          else
            C = 108584422*2^(-26)+120580432*2^(-53);
            Cinf = -131340486*2^(-81)-133428324*2^(-109);
            Csup = -131340486*2^(-81)-133428320*2^(-109);
          end
        case 'log2'
          % 0.693147180559945309417232121458176568075500134360 ...
          if nargout<=1
            Cinf = 93032639*2^(-27)+99906526*2^(-54);
            Csup = 93032639*2^(-27)+99906528*2^(-54);
          else
            C = 93032639*2^(-27)+99906526*2^(-54);
            Cinf = 112142222*2^(-82)+108200062*2^(-109);
            Csup = 112142222*2^(-82)+108200064*2^(-109);
          end
        case 'log10'
          % 2.3025850929940456840179914546843642076011014886287 ...
          if nargout<=1
            Cinf = 77261934*2^(-25)+124430890*2^(-52);
            Csup = 77261934*2^(-25)+124430892*2^(-52);
          else
            C = 77261934*2^(-25)+124430892*2^(-52);
            Cinf = -131214162*2^(-79)-87723944*2^(-107);
            Csup = -131214162*2^(-79)-87723940*2^(-107);
          end
        case 'sqrt2'
          % 1.414213562373095048801688724209698078569671875376948073176679 ...
          if nargout<=1
            Cinf = 94906265*2^(-26)+83785624*2^(-53);
            Csup = 94906265*2^(-26)+83785626*2^(-53);
          else
            C = 94906265*2^(-26)+83785626*2^(-53);
            Cinf = -116870404*2^(-80)-124045484*2^(-107);
            Csup = -116870404*2^(-80)-124045482*2^(-107);
          end
        case 'sqrt3'
          % 1.732050807568877293527446341505872366942805253810 ...
          if nargout<=1
            Cinf = 116235962*2^(-26)+92588704*2^(-56);
            Csup = 116235962*2^(-26)+92588720*2^(-56);
          else
            C = 116235962*2^(-26)+92588704*2^(-56);
            Cinf = 121316724*2^(-80)+95293200*2^(-109);
            Csup = 121316724*2^(-80)+95293208*2^(-109);
          end
        case 'sqrt5'
          % 2.2360679774997896964091736687312762354406183596115257242708 ...
          if nargout<=1
            Cinf = 75029990*2^(-25)+120580430*2^(-52);
            Csup = 75029990*2^(-25)+120580432*2^(-52);
          else
            C = 75029990*2^(-25)+120580432*2^(-52);
            Cinf = -131340486*2^(-80)-133428324*2^(-108);
            Csup = -131340486*2^(-80)-133428320*2^(-108);
          end
        otherwise
          if nargout==2
            [c,err] = str2intval(x);
          else
            c = str2intval(x);
          end
          setround(rndold)
          return
      end
      if nargout==2
        c = C;
        err.complex = 0;
        err.inf = Cinf;
        err.sup = Csup;
        err.mid = [];
        err.rad = [];
        err = class(err,'intval');
      else
        c.complex = 0;
        c.inf = Cinf;
        c.sup = Csup;
        c.mid = [];
        c.rad = [];
        c = class(c,'intval');
      end
      
    elseif isa(x,'polynom')     % 1 parameter, polynom
      c = polynom(intval(vector(x)));

    elseif isa(x,'gradient')    % 1 parameter, gradient
      if isa(x.x,'intval')
        c = x;
      else
        c = gradient(x,'gradientintval');
      end

    elseif isa(x,'hessian')    % 1 parameter, Hessian
      if isa(x.x,'intval')
        c = x;
      else
        c = hessian(x,'hessianintval');
      end
      
    elseif isa(x,'taylor')    % 1 parameter, Taylor
      if isa(x.t,'intval')
        c = x;
      else
        c = taylor(x,'taylorintval');
      end
      
    elseif isa(x,'slope')       % 1 parameter, slope
      c = x;                    % slopes are always of type intval

    else
      error('invalid call of intval constructor')
    end

  elseif nargin==3                  % internal type cast
    if isequal(type,'midrad')       % internal use only
      % input real or complex, rad is real nonnegative
      c.complex = 1;
      c.inf = [];
      c.sup = [];
      if prod(size(x))==1
        c.mid = x*ones(size(y));
      else
        c.mid = x;
      end
      % careful: just 0 is not sparse and may cause tremendous memory need
      if isequal(y,0)
        c.rad = 0;
      else
        if prod(size(y))==1
          if issparse(x)
            c.rad = sparse(repmat(y,size(x)));
          else
            c.rad = repmat(y,size(x));
          end
        else
          if ~isequal(size(c.mid),size(y))
            error('invalid call of intval constructor')
          else
            c.rad = y;
          end
        end
      end
    elseif isequal(type,'infsup')   % internal use only
      % input is real and inf<=sup
      if ~isreal(x) | ~isreal(y)
        disp('*************************************************************')
        disp('*************************************************************')
        disp('*************************************************************')
        disp('This should *never* happen!')
        disp('Please report circumstances to the author, rump@tuhh.de')
        disp('Thanks for your cooperation.')
        disp('*************************************************************')
        disp('*************************************************************')
        disp('*************************************************************')
        error('invalid use of intval(x,y,''infsup'')')
      end
      c.complex = 0;
      if issparse(x) | issparse(y)
        if prod(size(x))==1
          if isequal(x,0)
            [m,n] = size(y);
            x = sparse([],[],[],m,n);
          else
            % cures Matlab 6.5 bug (no problem since x is scalar)
            x = sparse(repmat(full(x),size(y)));
          end
        elseif prod(size(y))==1
          if isequal(y,0)
            [m,n] = size(x);
            y = sparse([],[],[],m,n);
          else
            y = repmat(y,size(x));
          end 
        end
      else
        if prod(size(x))==1
          x = repmat(x,size(y));
        elseif prod(size(y))==1
          y = repmat(y,size(x));
        end
      end
      c.inf = x;
      c.sup = y;
      c.mid = [];
      c.rad = [];
    else
      error('invalid call of intval constructor')
    end
    c = class(c,'intval');
  elseif nargin==0                % for save/load
    c.complex = 0;
    c.inf = [];
    c.sup = [];
    c.mid = [];
    c.rad = [];
    c = class(c,'intval');
  else
    error('invalid call of intval constructor')
  end

  % avoid Matlab 6.5f bug: 
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if isa(c,'intval')
    if c.complex
      if prod(size(c.mid))==1
        c.mid = full(c.mid);
        c.rad = full(c.rad);
      end
    else
      if prod(size(c.inf))==1
        c.inf = full(c.inf);
        c.sup = full(c.sup);
      end
    end
  end
    
  if rndold
    setround(rndold)
  end
