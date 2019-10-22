function r = gradient(a,str)
%GRADIENT     Gradient class constructor for interval input
%
%   r = gradient(a)
%
%An explicit call of the constructor is only necessary to initialize
%  a constant to be of type gradient. Otherwise, any operation
%  with a dependent variable produces a result of type gradient.
%

% gradient.x always stored as column vector with corresponding
%   matrix of gradients gradients.dx
% actual size of gradient is size(a.x)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump
% modified 12/18/02     S.M. Rump  Comment corrected
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    extra parameter str for 'random' (only for test purposes)

% modified 08/01/04     S.M. Rump  undocumented: accepts individual specification of variables
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_GRADIENT_NUMVAR = getappdata(0,'INTLAB_GRADIENT_NUMVAR');

  if INTLAB_GRADIENT_NUMVAR==0
    error('no dependent variables initialized for use of gradient')
  end

  dummy.init = a;
  if nargin==1
    r = gradient(dummy,'gradient');
  else
    r = gradient(dummy,str);
  end
