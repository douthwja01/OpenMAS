function u = diag(a,k)
%DIAG         Implements  diag(a,k)  for slopes
%
%   u = diag(a,k)
%
% functionality as Matlab function diag for matrices
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  if nargin==1
    k = 0;
  end

  if length(a.size)>2
    error('function slope/diag only for vectors and matrices')
  end

  if ( a.size(1)==1 ) | ( a.size(2)==1 )
    index = diag(ones(1,prod(a.size)),k);
    u.size = size(index);
    if issparse(a.r)
      u.r = intval(sparse([],[],[],prod(size(index)),INTLAB_SLOPE.NUMVAR+1,0));
      u.s = intval(sparse([],[],[],prod(size(index)),INTLAB_SLOPE.NUMVAR,0));
    else
      u.r = intval(zeros(prod(size(index)),INTLAB_SLOPE.NUMVAR+1));
      u.s = intval(zeros(prod(size(index)),INTLAB_SLOPE.NUMVAR));
    end
    u.r(index~=0,:) = a.r;
    u.s(index~=0,:) = a.s;
  else
    index = diag( reshape( 1:prod(a.size) , a.size ) , k );
    u.size = size(index);
    u.r = a.r( index(:) , : );
    u.s = a.s( index(:) , : );
  end

  u = class(u,'slope');
