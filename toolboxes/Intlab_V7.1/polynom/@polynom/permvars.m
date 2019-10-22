function p = permvars(p,perm)
%PERMVARS     Permute variables of p according to perm
%
%   q = permvars(p,perm)
%
%For given polynomial p in n variables, perm is a permutation of 1:n.
%On return, q is the polynomial p with variables permuted according to perm.
%
%   q = permvars(p,'lex')
%
%sorts variables of p in lexicographical order
%

% written  09/17/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  n = size(p.e,2);
  if n==1                           % nothing to do
    return
  end

  if isequal(perm,'lex')            % lexicographical order
    [p.v,index] = sort(p.v);
    p.e = p.e(:,index);
    [p.e,index] = sortrows(p.e);
    p.e = flipud(p.e);
    p.c = flipud(p.c(index));
    return
  end

  if ~isequal(sort(perm(:)),(1:n)')
    error('second argument is no permutation of 1:n, n=numvars(p)')
  end

  if n~=1
    p.e = p.e(:,perm);
    p.v = p.v(perm);
  end
