function res = longinit(param,see)
%LONGINIT     Initialization and defaults of long toolbox
%
%   longinit(param,see)
%
%possible values for param:
%
%  'Init'                    Initialize INTLAB constants (for startup file)
%
%Default setup of long arithmetic:
%  'WithErrorTerm'           Every long number with error term
%  'WithoutErrorTerm'        Pure long number arithmetic, no error terms
%  'ErrorTerm'               res = 'WithErrorTerm'
%                                  'WithoutErrorTerm'
%
%Long numbers with error term implement a midpoint/radius arithmetic for
%  long number with long midpoint and double radius. The radii are scaled
%  such that large exponents and (very) long numbers are possible.
%Default is with error term. Mantissas are shortened according to error term.
%
%When switching to 'WithErrorTerm', further calculations assume long number
%  to have an error term. To avoid errors use addlongerror.
%

% written  12/30/98     S.M. Rump
% modified 10/01/12     S.M. Rump  comment
%

%Large beta implies higher precision but smaller max. length pmax to be
%handled by one convolution. Maximum representable integer bitmax is
%2^53-1, therefore pmax=floor(bitmax/beta/(beta-1)):
%
%  logbeta              25     24     23     22     21     20     19
%  pmax                  8     32    128    512   2048   8192  32768
%  precision in bit    200    768   2944  11264  43008 163840 622592
%  precision in dec     60    231    886   3390  12946  49320 187418
%
%Lengths exceeding pmax slow down computation significantly
%
%A corresponding message is printed unless when changing a default mode
%  unless it is explicitly suppressed with the (optional) parameter see=0.
%

% written  10/16/98     S.M. Rump
% modified 07/28/02     S.M. Rump  parameter see added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%
  
  if nargin==1
    see = 1;
  end

switch lower(param)

  case 'init'

    INTLAB_LONG_LOGBETA = 23;      % long precision, <=25
    if INTLAB_LONG_LOGBETA > 25
      error('invalid initialization of INTLAB_LONG_LOGBETA in longinit');
    end
    INTLAB_LONG_BETA = pow2(INTLAB_LONG_LOGBETA);
    setappdata(0,'INTLAB_LONG_BETA',INTLAB_LONG_BETA);
    setappdata(0,'INTLAB_LONG_LOGBETA',INTLAB_LONG_LOGBETA);

    longprecision(0);              % set LongPrecision to minimum precision

    % default with error terms (midpoint/radius interval arithmetic
    setappdata(0,'INTLAB_LONG_ERROR',1);         

  case 'witherrorterm'
    setappdata(0,'INTLAB_LONG_ERROR',1);
    if see
      disp('===> Long arithmetic computations with error terms, i.e. valid long error bounds')
    end
    return

  case 'withouterrorterm'
    setappdata(0,'INTLAB_LONG_ERROR',0);
    if see
      disp('===> Long arithmetic computations without error terms')
    end
    return

  case 'errorterm'
    INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');
    if INTLAB_LONG_ERROR
      res = 'WithErrorTerm';
    else
      res = 'WithoutErrorTerm';
    end
    return

  otherwise
    error('longinit called with invalid argument')

  end
