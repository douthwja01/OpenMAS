function disp_(c)
%DISP_        Display of interval hessians in "_" notation
%
%   disp_(c)
%

% written  06/04/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  if nargin<2
    name = inputname(1);
    if isempty(name)                    % happens for display(hessianinit(random))
      name = 'ans';
    end
  end
  
  if isa(c.x,'intval')
    display_gen(c,name,'disp_')
  else
    display_gen(c,name,'display_')
  end

  