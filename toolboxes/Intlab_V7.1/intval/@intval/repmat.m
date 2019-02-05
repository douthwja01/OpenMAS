function a = repmat(a,varargin)
%REPMAT       Implements  repmat(a)  for interval data
%
%Functionality as in Matlab.
%

% written  02/17/08     S.M. Rump 
%

  if a.complex                  % complex input
    if isequal(a.rad,0)         % zero radius
      a = intval(repmat(a.mid,varargin{:}),0,'midrad');
    else                        % nonzero radius
      a = intval(repmat(a.mid,varargin{:}),repmat(a.rad,varargin{:}),'midrad');
    end
  else                          % real input
    a = intval(repmat(a.inf,varargin{:}),repmat(a.sup,varargin{:}),'infsup');
  end
  