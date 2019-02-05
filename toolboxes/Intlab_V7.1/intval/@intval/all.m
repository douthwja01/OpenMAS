function res = all(a,dim)
%ALL          Like Matlab function "all" for intval
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
    if a.complex                % complex input
      if isequal(a.rad,0)
        res = all(a.mid);
      else
        res = all(a.mid) | all(a.rad);
      end
    else                        % real input
      res = all(a.inf) | all(a.sup);
    end
  else
    if a.complex                % complex input
      if isequal(a.rad,0)
        res = all(a.mid,dim);
      else
        res = all(a.mid,dim) | all(a.rad,dim);
      end
    else                        % real input
      res = all(a.inf,dim) | all(a.sup,dim);
    end
  end
