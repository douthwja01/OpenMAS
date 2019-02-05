function res = sparsehessian(N,see)
%SPARSEHESSIAN   Number of unknowns from which hessians are stored sparse
%
%   N = sparsehessian         gets current value, default is 10
%   res = sparsehessian(N)    forces the gradient and hessian part 
%                               of hessians to be sparse for N and more unknowns
%
%The values N=0 and N=inf imply always or never to use sparse storage, respectively. 
%The default is N=10 (set in startintlab), i.e. to use sparse storage for ten and 
%  more variables. 
%
%When changing the current value, a corresponding message will be printed.
%To suppress the message, use
%
%   res = sparsehessian(N,0) 
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/05/05     S.M. Rump  comment corrected
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
  if ( nargin==0 )                  % getting current width
    res = INTLAB_HESSIAN_SPARSE;
  else                              % setting current width
    if ( round(N)~=N ) | ( N<0 )
      error('Invalid input for sparsehessian')
    end
    INTLAB_HESSIAN_SPARSE = N;
    res = INTLAB_HESSIAN_SPARSE;
    setappdata(0,'INTLAB_HESSIAN_SPARSE',INTLAB_HESSIAN_SPARSE);
    if nargin==1
      see = 1;
    end
    if see
      if INTLAB_HESSIAN_SPARSE==0
        disp('===> Hessian derivatives always stored sparse')
      elseif INTLAB_HESSIAN_SPARSE==inf
        disp('===> Hessian derivatives always stored full')
      else
        disp(['===> Hessian derivatives stored sparse for ' int2str(N) ' and more unknowns'])
      end
    end
  end
  