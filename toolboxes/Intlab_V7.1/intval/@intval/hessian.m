function r = hessian(a,str)
%HESSIAN      Hessian class constructor for interval input
%
%   r = hessian(a)
%
%An explicit call of the constructor is only necessary to initialize
%  a constant to be of type hessian. Otherwise, any operation
%  with a dependent variable produces a result of type hessian.
%

% written  12/18/02     S.M. Rump
% modified 03/06/04     S.M. Rump  Hessians allowed for non-scalars
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    extra parameter str for 'random' (only for test purposes)
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_HESSIAN_NUMVAR = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  if INTLAB_HESSIAN_NUMVAR==0
    error('no dependent variables initialized for use of hessian')
  end

  dummy.init = a;
  if nargin==1
    r = hessian(dummy,'hessian');
  else
    if ~isequal(str,'random')
      error('invalid call of @intval/hessian')
    end
    r = hessian(dummy,str);
  end
