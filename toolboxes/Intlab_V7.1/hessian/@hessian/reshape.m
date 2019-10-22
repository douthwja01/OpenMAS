function a = reshape(a,varargin)
%RESHAPE      Reshape for hessian vectors/matrices
%
%   r = reshape(a,vector)  or  r = reshape(a,n1,n2,...)
%
% functionality as Matlab function reshape
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a.x = reshape(a.x,varargin{:});
