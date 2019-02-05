function k = longprecision(k)
%LONGPRECISION  Sets/gets maximum precision for long number arithmetic
%
%   k = longprecision       gets precision in decimals
%   longprecision(k)        sets precision (approximately) in decimals
%                             and gives back actual precision
%   longprecision(0)        sets precision to minimum precision necesary
%                             to store doubles w/o error
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');

  if ( nargout==0 ) & ( nargin==1 )        % set precision

    RequiredPrecision = ceil( k*log(10)/log(2)/INTLAB_LONG_LOGBETA );
    MimimumPrecision = ceil(52/INTLAB_LONG_LOGBETA) + 1 ;
    INTLAB_LONG_PRECISION = max( RequiredPrecision , MimimumPrecision );
    k = floor( INTLAB_LONG_PRECISION*INTLAB_LONG_LOGBETA*log10(2) );
    setappdata(0,'INTLAB_LONG_PRECISION',INTLAB_LONG_PRECISION);

  elseif ( nargout<=1 ) & ( nargin==0 )    % get precision

    INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
    k = floor( INTLAB_LONG_PRECISION*INTLAB_LONG_LOGBETA*log10(2) );

  else
    error('invalid call of longprecision')
  end
