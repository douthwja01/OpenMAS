function a = sum(a,dim)
%SUM          Implements  sum(a,dim)  for Taylor
%
%   c = sum(a,dim)
%
% parameter dim optional, functionality as Matlab function sum
%

% written  05/21/09     S.M. Rump
%

  m = a.size(1);
  n = a.size(2);
  k = size(a.t,1);
  if nargin==1
    if m==1
      dim = 2;
    else
      dim = 1;
    end
  end

  if ( m==1 ) | ( n==1 )    % vector sum
    if  ~( ( ( m==1 ) & ( dim==1 ) ) | ( ( n==1 ) & ( dim==2 ) ) )
      a.t = sum(a.t,2);
    end
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if dim==1                             % sum of columns
    index = reshape(1:prod(a.size),m,n)';
    a.size = [1 n];
    a.t = reshape(sum(reshape(a.t(:,index),k*n,m),2),k,n);   % sum of rows of transposed
  else                                  % sum of rows
    a.size = [m 1];
    a.t = reshape(sum(reshape(a.t,k*m,n),2),k,m);
  end
  
  if rndold
    setround(rndold)
  end
