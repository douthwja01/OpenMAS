function a = reshape(a,varargin)
%RESHAPE      Reshape for Taylor vectors/matrices
%
%   r = reshape(a,vector)  or  r = reshape(a,n1,n2,...)
%
% functionality as Matlab function reshape
%

% written  05/21/09     S.M. Rump
% modified 09/28/12     S.M. Rump  performance improvement
%

  index = reshape(1:prod(a.size),varargin{:});
  a.size = size(index);
