function a = sum(a,dim)
%SUM          Implements  sum(a,dim)  for gradients
%
%   c = sum(a,dim)
%
% parameter dim optional, functionality as Matlab function sum
%

% written  10/16/98     S.M. Rump
% modified 11/03/03     S.M. Rump  performance improvement
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/05/05     S.M. Rump  improved performance
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
      a.dx = sum(a.dx);
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

  a.x = sum(a.x,dim);       % matrix sum
  if dim==1                 % column sum
    a.dx = sparse(repmat(1:n,m,1),1:m*n,1)*a.dx;
  else                      % row sum
    a.dx = repmat(speye(m),1,n)*a.dx;   % thanks to J. Kubitz
  end
  
  if rndold
    setround(rndold)
  end
