function i = end(p,k,n)
%END          Overloaded functions end, specifies last index
%
%This indexing is ambiguous for polynomials ( p(end)=p_n or p_0 ) and therefore forbidden

% written  10/03/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  error('indexing with end ambiguous and therefore prohibited for polynomials')
