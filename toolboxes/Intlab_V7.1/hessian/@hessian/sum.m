function a = sum(a,dim)
%SUM          Implements  sum(a,dim)  for hessians
%
%   c = sum(a,dim)
%
% parameter dim optional, functionality as Matlab function sum
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  [m n] = size(a.x);
  if nargin==1
    if m==1
      dim = 2;
    else
      dim = 1;
    end
  end

  if ( m==1 ) | ( n==1 )    % vector sum
    if  ~( ( ( m==1 ) & ( dim==1 ) ) | ( ( n==1 ) & ( dim==2 ) ) )
      a.x = sum(a.x);
      a.dx = sum(a.dx,2);
      a.hx = sum(a.hx,2);
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

  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  a.x = sum(a.x,dim);       % matrix sum
  if dim==1                 % column sum    
    index = reshape(1:m*n,m,n)';
    a.dx = reshape(sum(reshape(a.dx(:,index),n*N,m),2),N,n);
    a.hx = reshape(sum(reshape(a.hx(:,index),n*N*N,m),2),N^2,n);
  else                      % row sum
    a.dx = reshape(sum(reshape(a.dx,m*N,n),2),N,m);
    a.hx = reshape(sum(reshape(a.hx,m*N*N,n),2),N^2,m);
  end
  
  if rndold
    setround(rndold)
  end
