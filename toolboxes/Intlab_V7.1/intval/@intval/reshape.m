function a = reshape(a,varargin)
%RESHAPE      Reshape for interval vectors/matrices
%
%   r = reshape(a,vector)  or  r = reshape(a,n1,n2,...)
%
% functionality as Matlab function reshape
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 06/16/13     S.M. Rump  performance improvement
%

  if a.complex
    a.mid = reshape(a.mid,varargin{:});
    if ~isequal(a.rad,0)
      a.rad = reshape(a.rad,varargin{:});
    end
  else
    a.inf = reshape(a.inf,varargin{:});
    a.sup = reshape(a.sup,varargin{:});
  end
