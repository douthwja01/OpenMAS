function X = sparse(varargin)
%SPARSE       type cast to sparse interval matrix
%
%Call
%   Y = sparse(X)
%is simple type cast of X to sparse Y
%
%Call
%   Y = sparse(i,j,s,m,n,nzmax)
%has same functionality as Matlab/sparse for intval quantity s; 
%  produces sparseintval quantity Y.
%

% written  06/21/99     S.M. Rump
% modified 09/01/00     S.M. Rump  improved performance
% modified 07/08/02     S.M. Rump  zero radius
% modified 12/08/02     S.M. Rump  adapted to Matlab/sparse
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/03/12     S.M. Rump  sparse up to 2 dimensions
%

  if length(varargin)==1          % simple type cast
    X = varargin{1};
    if X.complex
      if length(size(X.mid))>2
        error('sparse arrays only up to 2 dimensions')
      end
      X.mid = sparse(X.mid);
      if ~isequal(X.rad,0)
        X.rad = sparse(X.rad);
      end
    else
      if length(size(X.inf))>2
        error('sparse arrays only up to 2 dimensions')
      end
      X.inf = sparse(X.inf);
      X.sup = sparse(X.sup);
    end
  else                            % Matlab functionality
    s = varargin{3};
    if s.complex
      X.complex = 1;
      X.inf = [];
      X.sup = [];
      varargin{3} = s.mid;
      X.mid = sparse(varargin{:});
      if isequal(s.rad,0)
        X.rad = 0;
      else
        varargin{3} = s.rad;
        X.rad = sparse(varargin{:});
      end
    else
      X.complex = 0;
      varargin{3} = s.inf;
      X.inf = sparse(varargin{:});
      varargin{3} = s.sup;
      X.sup = sparse(varargin{:});
      X.mid = [];
      X.rad = [];
    end
    X = class(X,'intval');
  end
  