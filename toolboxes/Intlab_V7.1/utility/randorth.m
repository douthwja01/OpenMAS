function res = RandOrth(n,m);
%RANDORTH     Random orthogonal matrix
%
%   res = randorth(n);
%

% written   1/21/95     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 11/29/10     S.M. Rump  fast version for larger dimension
% modified 05/28/11     S.M. Rump  random row and column permutations
%                                    and using Householder for large dimension
% modified 05/29/11     S.M. Rump  rectangular matrices
% modified 08/26/12     S.M. Rump  global variables removed
% modified 12/06/12     S.M. Rump  lower case
%

  e = 1e-30;
  if 1+e==1-e                   % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  if nargin==2
    k = max(m,n);
  else
    k = n;
  end

  if k<50                       % small k: generate new orthogonal matrix
    res = orth(randn(k));
  else                          % large k: generate new or permute old orth. matrix
    INTLAB_ORTH = getappdata(0,'INTLAB_ORTH');
    if ( length(INTLAB_ORTH)<k ) | isempty(INTLAB_ORTH{k})
      % generate new matrix of dimension k
      if k<500
        INTLAB_ORTH{k} = orth(randn(k));
        setappdata(0,'INTLAB_ORTH',INTLAB_ORTH);
      else
        A = eye(k);
        for i=1:50
          v = randn(k,1);
          A = A - ( (2/(v'*v) )*v ) * (v'*A);
        end
        INTLAB_ORTH{k} = A;
        setappdata(0,'INTLAB_ORTH',INTLAB_ORTH);
      end
      res = INTLAB_ORTH{k};
    else
      res = INTLAB_ORTH{k}(randperm(k),randperm(k));
    end
  end
  if nargin==2
    if n>m
      res = res(:,1:m);
    elseif m>n
      res = res(1:n,:);
    end
  end
  
  if rndold
    setround(rndold)
  end
