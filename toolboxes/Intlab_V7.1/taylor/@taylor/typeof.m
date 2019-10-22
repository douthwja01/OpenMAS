function t = typeof(a)
%TYPEOF       Type of a
%
%   t = typeof(a)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  05/22/09     S.M. Rump  taylor added
%

% taylor\@taylor\typeof:  a  must be taylor
  if isa(a.t,'intval')
    t = 'taylorintval';
  else
    t = 'taylor';
  end
