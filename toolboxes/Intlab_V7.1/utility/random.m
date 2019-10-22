function r = random(varargin)
%RANDOM       Random numbers in [min,max],  default [-1,+1]
%
% Calling conventions as function  rand
%   random         single random number
%   random(n)      n x n random matrix
%   random(m,n)    m x n random matrix           same as random([m,n])
%   random(m,n,k)  m x n x k random matrix       same as random([m,n,k])
%             .........
%
% If the last argument is a cell array {min,max} of length 2,
%   then random numbers uniformly within [min,max] are generated
% default is [-1,1]
%
%   random(3,{.1,.2})   3 x 3 random matrix, entries between .1 and .2
%

% written  02/08/91     S.M. Rump
% modified 03/18/98     S.M. Rump    user defined interval
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/19/08     S.M. Rump  point range
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  n = length(varargin);
  if n>0 & ~iscell(varargin{1}) & isnan(varargin{1})      % call from randomc
    if n>1
      varargin = varargin{2};
      n = length(varargin);
    else
      n = 0;
    end
  end

  if n>0
    if iscell(varargin{n})
      min = varargin{n}{1};
      max = varargin{n}{2};
      if min>max
        error('invalid margins for random numbers')
      end
      n = n-1;
    else
      min = -1;
      max = +1;
    end
  else
    min = -1;
    max = +1;
  end

  if n>0
    if length(varargin{1})>1
      s = [ '([' num2str(varargin{1}) '])' ];
    else
      s = '(';
      for i=1:n
        s = [s num2str(varargin{i}) ','];
      end
      s = [s(1:end-1) ')'];
    end
  else
    s=[];
  end

  eval(['r = rand' s ';']);
  r = r * ( max - min ) + min ;
  
  if rndold
    setround(rndold)
  end
