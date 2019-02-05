function C = ProdKL(A,B,K,L)
%PRODKL       Matrix product in approximately K-fold precision stored in L results
%
%   C = ProdKL(A,B,K,L)
%
% Input A or B or both may be cell arrays, 1 <= L <= K. 
% Output C is cell array iff L>1, default is L=1.
% Simple application of SumKL. All input must be real.
% Relative error of the result is of order  eps^L + eps^K cond(product) , see
%   S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
%      Japan J. Indust. Appl. Math. (JJIAM), 26:249-277, 2009.
%
%Reference implementation! Slow due to interpretation!
%

% written  02/17/08     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, default L=1
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==3
    L = 1;
  end
  
  if iscell(A)              % input A cell array
    n = size(A{1},1);
    lenA = length(A);
    AA = zeros(n,lenA*n);
    for i=1:lenA
      AA(:,(i-1)*n+1:i*n) = A{i};
    end
    if iscell(B)            % both A and B cell arrays
      lenB = length(B);
      m = lenA*n;
      AAA = zeros(n,lenB*m);
      BBB = zeros(lenB*m,n);
      for i=1:lenB
        AAA(:,(i-1)*m+1:i*m) = AA;
        BBB((i-1)*m+1:i*m,:) = repmat(B{i},lenA,1);
      end
      A = AAA;
      B = BBB;
    else
      A = AA;
      B = repmat(B,lenA,1);
    end
  else
    n = size(A,1);
  end
  if iscell(B)              % input B cell array
    lenB = length(B);
    A = repmat(A,1,lenB);
    BB = zeros(lenB*n,n);
    for i=1:lenB
      BB((i-1)*n+1:i*n,:) = B{i};
    end
    B = BB;
  end
  if L==1
    C = zeros(n);
  else
    for i=1:L
      C{i} = zeros(n);
    end
  end
  for i=1:n
    for j=1:n
      [x,y] = TwoProduct(A(i,:),B(:,j)');
      res = SumKL([x y],K,L);
      if L==1
        C(i,j) = res;
      else
        for k=1:L
          C{k}(i,j) = res(k);
        end
      end
    end
  end
  
  if rndold
    setround(rndold)
  end
