function Err = errornormalize(Err)
%ERRORNORMALIZE  Normalization of error term
%
%For internal use only

%
%  Err = errornormalize(Err)
%
%with Err = Err.mant * beta^Err.exp and Err.mant normalized to
%approximately 1.
%

% written  12/30/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');

  % adapt error mantissa and exponent such that Err.mant around 1
  index = ( Err.mant~=0 );
  if any(index)
    q = floor( log( Err.mant(index) ) / log( INTLAB_LONG_BETA ) );
    Err.mant(index) = Err.mant(index) ./ INTLAB_LONG_BETA.^q;
    Err.exp(index) = Err.exp(index) + q;
  end
  Err.exp(~index) = 0;
