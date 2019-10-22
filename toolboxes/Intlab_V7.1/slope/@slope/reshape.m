function a = reshape(a,varargin)
%RESHAPE      Reshape for slope vectors/matrices
%
%   r = reshape(a,vector)  or  r = reshape(a,n1,n2,...)
%
% functionality as Matlab function reshape
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/12     S.M. Rump  performance improvement
%

  index = reshape(1:prod(a.size),varargin{:});
  a.size = size(index);
