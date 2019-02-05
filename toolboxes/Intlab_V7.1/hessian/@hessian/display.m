function display(c,name)
%DISPLAY      Command window display of hessian
%
%First and second derivative of Hessians are stored sparse if the number of dependent variables
%exceeds a certain constant (or details, see hessianinit).
%
%To change the display of a hessian variable "u" use
%
%   full(u)   or   sparse(u) .
%
%Careful, if the number of independent variables is large and full storage of second derivatives is
%chosen, output may become large.
%

%Second parameter name for internal purposes
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 08/26/12     S.M. Rump  global variables removed
%

  if nargin<2
    name = inputname(1);
    if isempty(name)                    % happens for display(hessianinit(random))
      name = 'ans';
    end
  end
  
  display_gen(c,name,'display_')
  

