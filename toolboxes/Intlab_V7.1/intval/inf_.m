function r = inf_(a)
%INF_         Implements  inf(a)  (cures problems with inf)
%
%   c = inf(a)
%
% On return, inf(a) <= alpha for all entries alpha in a
%
%************************************************************************
%********  due to conflict with internal variable inf (infinity)  *******
%********                    use function inf_                    *******
%************************************************************************
%

% Implementation for point arguments
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = a;
