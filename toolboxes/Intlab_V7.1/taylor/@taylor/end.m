function i = end(x,k,n)
%END          Overloaded functions end, specifies last index
%

% written  05/21/09     S.M. Rump
% modified 07/07/10     S.M. Rump  new design
% modified 08/26/12     S.M. Rump  global variables removed
%

  % This function is called by T(4:end) or T{2:end-1} or alike for a Taylor
  % variable T. Then meaning of "end" is different when called by () or {}.
  % However, this function "end" is called before entering the methods
  % subsref or subsasgn. I see no possibility to discriminate between these
  % two calls. For example, the calls T(2:end) and T(end-2:4) for an array
  % of length 4 yield the same index vector 2:4 passed by s.subs, but there
  % seems no way to recover this in subsref. Therefore, using "end" in an
  % index expression is only allowed for T(), not for T{}.
  
  % Should be  i = size(x.t,1)-1;  for reference by {}
  i = x.size(k);                % only correct for reference by ()
  setappdata(0,'INTLAB_TAYLOR_END',1);
