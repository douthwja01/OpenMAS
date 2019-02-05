function realimag(c)
%REALIMAG     Display real and imaginary part of interval gradients separately
%
%   realimag(c)
%

% written  10/16/98     S.M. Rump
% modified 12/18/02     S.M. Rump  imag corrected
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    name corrected
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  name = inputname(1);
  if isempty(name)                    % happens for display(gradientinit(random))
    name = 'ans';
  end
  
  if isreal(c.x)
    display(c,name)
  else
    if loose, disp(' '); end
    display(real(c),['real(' name ')'])
    if loose, disp(' '); end

    if loose, disp(' '); end
    display(imag(c),['imag(' name ')'])
    if loose, disp(' '); end
  end

  