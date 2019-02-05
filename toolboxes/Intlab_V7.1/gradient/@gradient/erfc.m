function a = erfc(a)
%ERFC         Gradient complementary error function erfc(a)
%

% written  05/31/13     S.M. Rump
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  N = getappdata(0,'INTLAB_GRADIENT_NUMVAR');
  % factorLB <= 2/sqrt(pi) <= factorUB
  INTLAB_STDFCTS_ERF = getappdata(0,'INTLAB_STDFCTS_ERF');
  factorLB = INTLAB_STDFCTS_ERF.TWO_SQRTPIINF;  % round to nearest
  factorUB = INTLAB_STDFCTS_ERF.TWO_SQRTPISUP;  % ~ 1.12

  % use full(a.x(:)): cures Matlab V6.0 bug
  % a=7; i=[1 1]; x=a(i), b=sparse(a); y=b(i)  yields row vector x but column vector y
  % ax is full anyway
  a.x = full(a.x);
  ax = exp(-a.x(:).^2);
  a.x = erfc(a.x);
  if issparse(a.dx)
    sizeax = size(a.dx,1);
    [ia,ja,sa] = find(a.dx);
    if isa(a.x,'intval')
      adx = times(ax(ia),sa,0);
      a.dx = intval(-factorUB,-factorLB,'infsup') * ...
          intval( sparse(ia,ja,adx.inf,sizeax,N) , sparse(ia,ja,adx.sup,sizeax,N) , 'infsup' );
    else
      a.dx = (-factorLB) * sparse(ia,ja,ax(ia).*sa,sizeax,N);
    end
  else
    if isa(a.x,'intval')
      ax = intval(-factorUB,-factorLB,'infsup') * ax;
      a.dx = a.dx .* ax(:,ones(1,N));
    else
      ax = (-factorLB) * ax;
      a.dx = a.dx .* ax(:,ones(1,N));
    end
  end
  
  if rndold
    setround(rndold)
  end
