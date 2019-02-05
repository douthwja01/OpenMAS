function rhx = adx2rhx(N,m,adx,bdx)
%For two N x m matrices adx,bdx fast computation of
%  rhx = sparse([],[],[],N^2,m)
%  for i=1:m
%    rhx(:,i) = reshape ( adx(:,i) * bdx(:,i).' , N^2,1 )
%  end
%
%If argument bdx is omitted, then bdx := adx .
%

  if nargin==4              % second factor bdx given

    ea = sum(spones(adx),1);
    eb = sum(spones(bdx),1);
    e = ea.*eb;
    m_ = m;
    
    % make sue no zero columns, makes index battle much easier
    if any(e==0)        
      Index = find(e);
      if isempty(Index)     % nothing to do
        rhx = sparse([],[],[],N^2,m);
        if isa(adx,'intval') | isa(bdx,'intval')
          rhx = intval(rhx);
        end
        return
      end
      adx = adx(:,Index);
      bdx = bdx(:,Index);
      ea = sum(spones(adx),1);
      eb = sum(spones(bdx),1);
      e = ea.*eb;
      m = length(Index);
      flag = 1;
    else
      flag = 0;
    end
    
    sa = nonzeros(adx);
    [ib,jb,sb] = find(bdx);
    
    % repeat eb(i) times each i in 1:m
    jr = 1:m;
    J = zeros(1,sum(eb));
    J(1) = 1;
    J(cumsum(eb(1:end-1))+1) = 1;
    J = cumsum(J);
    indexa = jr(J);
    % form index vector for first factor sadx
    [indexa,dummy,sadx] = find(adx(:,indexa));
    
    % repeat eb(i) times each i in ea
    ix = ea(J);
    % repeat ix(i) times each i in 1:length(sb)
    J = zeros(1,sum(ix));
    J(1) = 1;
    J(cumsum(ix(1:end-1))+1) = 1;
    J = cumsum(J);
    indexb = 1:length(sb);
    sbdx = reshape(sb(indexb(J)),size(sadx));
    
    % factors for outer product are sadx, sbdx
    % compute row vector for final result
    ir = indexa + reshape((ib(J)-1)*N,size(indexa));
    
  else            % second factor bdx equal to adx by default
    
    ea = sum(spones(adx),1);
    m_ = m;
    
    % make sue no zero columns, makes index battle much easier
    if any(ea==0)        
      Index = find(ea);
      if isempty(Index)     % nothing to do
        rhx = sparse([],[],[],N^2,m);
        if isa(adx,'intval')
          rhx = intval(rhx);
        end
        return
      end
      adx = adx(:,Index);
      ea = sum(spones(adx),1);
      m = length(Index);
      flag = 1;
    else
      flag = 0;
    end
    
    e = ea.*ea;
    [ia,ja,sa] = find(adx);
    
    % repeat ea(i) times each i in 1:m
    jr = 1:m;
    J = zeros(1,sum(ea));
    J(1) = 1;
    J(cumsum(ea(1:end-1))+1) = 1;
    J = cumsum(J);
    indexa = jr(J);
    % form index vector for first factor sadx
    [indexa,dummy,sadx] = find(adx(:,indexa));
    
    % repeat ea(i) times each i in ea
    ix = ea(J);
    % repeat ix(i) times each i in 1:length(sa)
    J = zeros(1,sum(ix));
    J(1) = 1;
    J(cumsum(ix(1:end-1))+1) = 1;
    J = cumsum(J);
    indexb = 1:length(sa);
    sbdx = reshape(sa(indexb(J)),size(sadx));
    
    % factors for outer product are sadx, sbdx
    % compute row vector for final result
    ir = indexa + reshape((ia(J)-1)*N,size(indexa));
    
  end
  
  % compute col vector for final result
  % repeat e(i) times each i in 1:m
  J = zeros(1,sum(e));
  J(1) = 1;
  J(cumsum(e(1:end-1))+1) = 1;
  jr = jr(cumsum(J));
  if flag
    jr = Index(jr);
  end
  
  % form final result
  if isa(sadx,'intval') | isa(sbdx,'intval')
    rhx = times(sadx,sbdx,0);
    if rhx.complex
      rhx = intval( sparse(ir,jr,rhx.mid,N^2,m_) , sparse(ir,jr,rhx.rad,N^2,m_) , 'midrad' );
    else
      rhx = intval( sparse(ir,jr,rhx.inf,N^2,m_) , sparse(ir,jr,rhx.sup,N^2,m_) , 'infsup' );
    end
  else
    rhx = sparse(ir,jr,sadx.*sbdx,N^2,m_);
  end
  
