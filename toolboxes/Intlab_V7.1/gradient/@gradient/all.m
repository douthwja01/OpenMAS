function res = all(a,dim)
%ALL          Like Matlab function "all" for gradient
%
%Call
%
%   L = all(A)
%   L = all(A,dim)
%
%Same functionality as Matlab/all for intval quantity A
%

% written  12/06/05     S.M. Rump
%

  if nargin==1
    res = all(a.x) | all( reshape(all(a.dx,2),size(a.x)) );
  else
    res = all(a.x,dim) | all( reshape(all(a.dx,2),size(a.x)) , dim );
  end
