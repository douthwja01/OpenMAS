function res = any(a,dim)
%ANY          Like Matlab function "any" for intval
%
%Call
%
%   L = any(A)
%   L = any(A,dim)
%
%Same functionality as Matlab/any for intval quantity A
%

% written  12/06/05     S.M. Rump
%

  if nargin==1
    if a.complex                % complex input
      if isequal(a.rad,0)
        res = any(a.mid);
      else
        res = any(a.mid) | any(a.rad);
      end
    else                        % real input
      res = any(a.inf) | any(a.sup);
    end
  else
    if a.complex                % complex input
      if isequal(a.rad,0)
        res = any(a.mid,dim);
      else
        res = any(a.mid,dim) | any(a.rad,dim);
      end
    else                        % real input
      res = any(a.inf,dim) | any(a.sup,dim);
    end
  end
