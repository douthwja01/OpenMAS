function n = degree(p,k);
%DEGREE       Degree of polynomial
%
%   n = degree(p);
%
%result:   degree of p for univariate polynomial
%          vector of degrees of p for multivariate polynomial
%
%   n = degree(p,k);
%
%result:   degree of p in the k-th variable; if p has less than k variables, then n=0.         
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin==2
    if k>size(p.e,2)
      n = 0;
    else
      n = max(p.e(:,k));
    end
  else
    if size(p.e,2)<=1       % univariate polynomial
      n = length(p.c)-1;
    else                    % multivariate polynomial
      n = max(p.e,[],1);
    end
  end
