function intlabtiming
% Timing of some INTLAB routines. Careful, may take very long due to large dimensions.
%

% written  10/14/12     S.M. Rump
%

TT = 10;                % Minimum time spend per test case to improve accuracy of timing
rr = midrad(1,1e-8);    % radius for interval data
disp(' ')
disp(['Vaio Intel i7 M640 2.8 GHz, Matlab ' version])

MatMul = 1;              % timing matrix multiplication
DenseLinSys = 1;         % timing dense linear systems
SparseLinSys = 1;        % timing sparse linear systems
OptNonlinSys = 1;        % timing nonlinear systems and optimization
IllcoLinSys = 1;         % timing extremely ill-conditoned linear systems (requires symbolic toolbox)
TimingApproxDot_ = 1;    % timing accurate dot product
TimingInclDot_ = 1;      % timing accurate dot product
AccuracyDot_ = 1;        % accuracy of accurate dot product
TimingLssresApprox1 = 1; % timing of approximate lssresidual, dense matrix
TimingLssresApprox2 = 1; % timing of approximate lssresidual, sparse matrix
TimingLssresApprox3 = 1; % timing of approximate lssresidual, matrix residual
TimingLssresIncl1 = 1;   % timing of verified lssresidual, dense matrix
TimingLssresIncl2 = 1;   % timing of verified lssresidual, sparse matrix
TimingLssresIncl3 = 1;   % timing of verified lssresidual, matrix residual
TimingDenseMatmul = 1;   % timing of verified lssresidual, matrix residual
TimingSparseMatmul = 1;  % timing of verified lssresidual, matrix residual
TimingVerifylss = 1;     % timing of verified lssresidual, matrix residual

if MatMul
  
  disp(' ')
  disp('Matrix multiplication, time [sec]')
  disp(' ')
  disp('               pure     verified    verified    verified')
  disp('  dimension   fl-pt        A*A       A*intA    intA*intA')
  disp('----------------------------------------------------------')

  N = [1000 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 4;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          switch m
            case 1, tic, A*A; t=t+toc;
            case 2, tic, intval(A)*A; t=t+toc;
            case 3, B = A*rr; tic, B*A; t=t+toc;
            case 4, A = A*rr; tic, A*A; t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff '  ' ff(T(m)) ];
        end
      end
    end
    
    disp(sprintf(['%8d ' fff],n,T))
    
  end
end

if DenseLinSys
  
  disp(' ')
  disp('Dense linear systems (up to factor 2 faster than previous version), time [sec]')
  disp(' ')
  disp('                  pure      verified      verified         verified')
  disp('  dimension      fl-pt         A\b      high acc. A\b     intA\intb')
  disp('---------------------------------------------------------------------')

  N = [1000 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 4;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          b = randn(n,1);
          switch m
            case 1, tic, A\b; t=t+toc;
            case 2, intvalinit('ImprovedResidual',0); tic, verifylss(A,b); t=t+toc;
            case 3, intvalinit('QuadrupleResidual',0); tic, verifylss(A,b); t=t+toc;
            case 4, intA = A*rr; intb = b*rr; tic, verifylss(intA,intb); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff '     ' ff(T(m)) ];
        end
      end
    end
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
end

if SparseLinSys
  
  disp(' ')
  disp('Sparse s.p.d. linear systems (approximately 10 nonzero elements per row), time [sec]')
  disp(' ')
  disp('               pure       symamd    verified    verified')
  disp('  dimension   fl-pt        fl-pt       A\b     intA\intb')
  disp('----------------------------------------------------------')

  N = [1000 2000 5000 10000 20000 50000];
  nmaxpureflpt = 10000;
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 4;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    intvalinit('ImprovedResidual',0);
    
    if n<=nmaxpureflpt
      Mmin = 1;             % all tests
    else
      Mmin = 2;             % avoid pure floating-point without preordering
    end
    for m=Mmin:M            % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = sprandn(n,n,3/n); 
          A = speye(n) + A'*A;
          v = symamd(A);
          B = A(v,v);
          b = randn(n,1);
          switch m
            case 1, tic, A\b; t=t+toc;
            case 2, tic, v = symamd(A); A(v,v)\b; t=t+toc;
            case 3, tic, verifylss(A,b); t=t+toc;
            case 4, intA = A*rr; intb = b*rr; tic, verifylss(intA,intb); t=t+toc;
          end
        end
        if phase==0                 % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                        % testing
          T(m) = t/K;
          fff = [ fff '  ' ff(T(m)) ];
        end
      end
    end
    
    if n<=nmaxpureflpt
      disp(sprintf(['%8d' fff],n,T))
    else
      disp(sprintf(['%8d        -   ' fff],n,T(2:end)))
    end
    
  end
end

if OptNonlinSys
  
  disp(' ')
  disp('Optimization and nonlinear systems, time [sec]')
  disp(' ')
  disp('             fminsearch    local   nonlin.system     local     maximum   verification ')
  disp('  dimension  fl-pt [sec]  minimum  verified [sec]   minimum   rel.error     pos.def.  ')
  disp('----------------------------------------------------------------------------------------')
  
  N = [50 100 300 1000 3000 10000];
  nmaxflpt = 1000;
  param = optimset('Display','off');
  for n=N
    approx = ones(n,1);
    if n<=nmaxflpt
      tic
      x = fminsearch(@(x) test_h(x,2),approx,param);
      t = toc;
      y = test_h(x,2);
    end
    tic
    X = verifynlss(@test_h,approx,'hSparseSPD',0,2);
    T1 = toc;    
    maxrelerr = max(relerr(X));
    tic
    Y = test_h(hessianinit(X),2);
    minimum = isspd(Y.hx);
    T2 = toc;
    if n<=nmaxflpt
      disp(sprintf(['%8d     ' ff(t) '   ' ff(y) ff(T1) '      ' ff(sup(Y.x)) '%9.1e    ' ff(T2)], ...
        n,t,y,T1,sup(Y.x),maxrelerr,T2))
    else
      disp(sprintf(['%8d          -          -      ' ff(T1) '      ' ff(sup(Y.x)) '%9.1e  ' ff(T2)], ...
        n,T1,sup(Y.x),maxrelerr,T2))
    end
  end
   
end


if IllcoLinSys
  
  disp(' ')
  disp('Ill-coonditioned linear systems, time [sec]')
  disp(' ')
  disp('                true       linear system    maximum  ')
  disp('  dimension  cond. number  verified [sec]  rel.error ')
  disp('-------------------------------------------------------')
  
  N = [10 20 30 40 50 100];
  for n=N
    A = invhilb(n);
    b = ones(n,1);
    cnd = TrueCond(A);
    tic
    X = verifylss(A,b,'illco');
    T = toc;    
    maxrelerr = max(relerr(X));
    disp(sprintf(['%8d     %10.1e     ' ff(T) '     %9.1e'],n,cnd,T,maxrelerr));
  end
   
end


if TimingApproxDot_

  disp(' ')
  disp('Accurate approximate dot product, time [sec]')
  disp(' ')
  disp('  dimension   t1       t2        t3    t2/t1     t3/t1    ')
  disp('------------------------------------------------------------')
  
  N = [500 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 3;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, b-A*x; t=t+toc;
            case 2, tic, lssresidual(A,x,b); t=t+toc;
            case 3, tic, Dot_(A,x,-1,b); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    for m=4:5
      T(m) = T(m-2)/T(1);
      fff = [ fff ff(T(m)) ];
    end
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if TimingInclDot_

  disp(' ')
  disp('Accurate verified dot product, time [sec]')
  disp(' ')
  disp('  dimension   t1       t2        t3      t2/t1     t3/t1    ')
  disp('--------------------------------------------------------------')
  
  N = [500 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 3;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, b-A*intval(x); t=t+toc;
            case 2, tic, lssresidual(A,x,intval(b)); t=t+toc;
            case 3, tic, Dot_(A,x,-1,b,-2); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    for m=4:5
      T(m) = T(m-2)/T(1);
      fff = [ fff ff(T(m)) ];
    end
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if AccuracyDot_

  disp(' ')
  disp('Accuracy dot product')
  disp('condition  median maximum    median maximum    median maximum')
  disp('---------------------------------------------------------------')
  
  n = 100;
  K = 100;
  for cnd=10.^[2 5 10 14]
    y1 = zeros(n,K);
    y2 = y1;
    y3 = y1;
    for k=1:K
      A = randmat(n,cnd);
      b = randn(n,1);
      x = A\b;
      y1(:,k) = rad(b-A*intval(x));
      y2(:,k) = rad(lssresidual(A,x,intval(b)));
      y3(:,k) = rad(Dot_(A,x,-1,b,-2));
    end
    y1 = y1(:);
    y2 = y2(:);
    y3 = y3(:);
    disp(sprintf('%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e', ...
      cnd,median(y1),max(y1),median(y2),max(y2),median(y3),max(y3)))
  end
   
end


if TimingLssresApprox1

  disp(' ')
  disp('Poor men''s approximate residual (lssresidual), dense matrix, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1 ')
  disp('----------------------------------------')
  
  N = [500 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, b-A*x; t=t+toc;
            case 2, tic, lssresidual(A,x,b); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if TimingLssresApprox2

  disp(' ')
  disp('Poor men''s approximate residual (lssresidual), sparse matrix, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1 ')
  disp('----------------------------------------')
  
  N = [5000 20000 50000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = sprand(n,n,20/n);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, b-A*x; t=t+toc;
            case 2, tic, lssresidual(A,x,b); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if TimingLssresApprox3

  disp(' ')
  disp('Poor men''s approximate residual (lssresidual), matrix residual, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1 ')
  disp('----------------------------------------')
  
  N = [500 2000 5000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = randn(n);
          R = inv(A);
          switch m
            case 1, tic, speye(n)-R*A; t=t+toc;
            case 2, tic, lssresidual(R,A,speye(n)); t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end

if TimingLssresIncl1

  disp(' ')
  disp('Poor men''s verified residual (lssresidual), dense matrix, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1      r1         r2 ')
  disp('-----------------------------------------------------------')
  
  N = [500 2000 5000];
  cnd = 1e8;
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    rr = zeros(1,M);
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        if phase==1
          r = zeros(1,K);
        end
        for k=1:K
          A = randmat(n,cnd);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, res = b-A*intval(x); t=t+toc;
            case 2, tic, res = lssresidual(A,x,intval(b)); t=t+toc;
          end
          r(k) = median(rad(res))/norm(A,inf);
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
      rr(m) = median(r);
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff '%10.1e %10.1e'],n,T,rr))
    
  end
   
end


if TimingLssresIncl2

  disp(' ')
  disp('Poor men''s verified residual (lssresidual), sparse matrix, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1      r1         r2 ')
  disp('-----------------------------------------------------------')
  
  N = [5000 20000 50000];
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    rr = zeros(1,M);
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        if phase==1
          r = zeros(1,K);
        end
        for k=1:K
          A = sprand(n,n,20/n);
          x = randn(n,1);
          b = A*x;
          switch m
            case 1, tic, res = b-A*intval(x); t=t+toc;
            case 2, tic, res = lssresidual(A,x,intval(b)); t=t+toc;
          end
          r(k) = median(rad(res))/norm(A,inf);;
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
      rr(m) = median(r);
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff '%10.1e %10.1e'],n,T,rr))
    
  end
   
end


if TimingLssresIncl3

  disp(' ')
  disp('Poor men''s verified residual (lssresidual), matrix residual, time [sec]')
  disp(' ')
  disp('  dimension   t1        t2      t2/t1      r1         r2 ')
  disp('-----------------------------------------------------------')
  
  N = [500 2000 5000];  
  cnd = 1e8;
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    rr = zeros(1,M);
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        if phase==1
          r = zeros(1,K);
        end
        for k=1:K
          A = randmat(n,cnd);
          R = inv(A);
          switch m
            case 1, tic, res = speye(n)-R*intval(A); t=t+toc;
            case 2, tic, res = lssresidual(R,A,intval(speye(n))); t=t+toc;
          end
          r(k) = median(rad(res(:)))/norm(A,inf);;
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
      rr(m) = median(r);
    end
    T(3) = T(2)/T(1);
    fff = [ fff ff(T(m)) ];
    
    disp(sprintf(['%8d' fff '%10.1e %10.1e'],n,T,rr))
    
  end
   
end


if TimingDenseMatmul

  disp(' ')
  disp('Dense matrix multiplication, time [sec]')
  disp(' ')
  disp('  dimension   t1       t2   ')
  disp('------------------------------')
  
  N = [500 2000 5000];  
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = midrad(randn(n),rand(n));
          switch m
            case 1, intvalinit('FastIVmult',0); tic, A*A; t=t+toc;
            case 2, intvalinit('SharpIVmult',0); tic, A*A; t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if TimingSparseMatmul

  disp(' ')
  disp('Sparse matrix multiplication, time [sec]')
  disp(' ')
  disp('  dimension   t1       t2   ')
  disp('------------------------------')
  
  N = [500 2000 5000];  
  for i=1:length(N)
    
    n = N(i);               % dimension
    M = 2;                  % number of test cases
    T = zeros(1,M);         % time for M test cases
    fff = '';               % format for output
    
    for m=1:M               % M test cases
      K = 1;
      for phase=0:1         % phase 0 timing, phase 1 testing
        t = 0;
        for k=1:K
          A = sprandn(n,n,0.05);
          [I,J,S] = find(A);
          A = A .* midrad(sparse(1),sparse(I,J,2*rand(size(S)).*abs(S)));
          switch m
            case 1, intvalinit('FastIVmult',0); tic, A*A; t=t+toc;
            case 2, intvalinit('SharpIVmult',0); tic, A*A; t=t+toc;
          end
        end
        if phase==0         % timing
          K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
        else                % testing
          T(m) = t/K;
          fff = [ fff ff(T(m)) ];
        end
      end
    end
    
    disp(sprintf(['%8d' fff],n,T))
    
  end
   
end


if TimingVerifylss

  N = [500 2000 5000];  
  disp(' ')
  disp('Dense linear system solver, time [sec]')
  disp(' ')
  for cnd=[1e8 1e14]
    disp(sprintf('cnd = %10.1e',cnd))
    disp('  dimension  t1        t2        t3      t2/t1     t3/t1      r1         r2         r3')
    disp('-------------------------------------------------------------------------------------------')
    for i=1:length(N)
      
      n = N(i);               % dimension
      M = 3;                  % number of test cases
      T = zeros(1,M);         % time for M test cases
      fff = '';               % format for output
      rr = zeros(1,M);
      
      for m=1:M               % M test cases
        K = 1;
        for phase=0:1         % phase 0 timing, phase 1 testing
          t = 0;
          if phase==1
            r = zeros(1,K);
          end
          for k=1:K
            A = randmat(n,cnd);
            b = randn(n,1);            
            switch m
              case 1, intvalinit('doubleresidual',0); tic, res = verifylss(A,b); t=t+toc;
              case 2, intvalinit('improvedresidual',0); tic, res = verifylss(A,b); t=t+toc;
              case 3, intvalinit('quadrupleresidual',0); tic, res = verifylss(A,b); t=t+toc;
            end
            r(k) = median(relerr(res(:)));
          end
          if phase==0         % timing
            K = min(100,max(1,ceil(TT/t)));  % spend at least TT time per test case
          else                % testing
            T(m) = t/K;
            fff = [ fff ff(T(m)) ];
          end
        end
        rr(m) = median(r);
      end
      for m=4:5
        T(m) = T(m-2)/T(1);
        fff = [ fff ff(T(m)) ];
      end
      
      disp(sprintf(['%8d' fff '%10.1e %10.1e %10.1e'],n,T,rr))
      
    end
  end
   
end


function str = ff(t)
% intelligent format: total length 10 characters, intelligent within %10.4f
% 4 digits following decimal point; factor 0.8 smoothens data near power of 10
L = max( 0 , min( 4 , 1 - floor(log10(0.8*t)) ));
if L==0
  D = 5;
else
  D = 6;
end
str = [ '%' int2str(D+L) '.' int2str(L) 'f' blanks(10-D-L) ];

function cnd = TrueCond(A)
% 100-Digits 2-norm condition number using vpa from symbolic toolbox
  Digits = 100;
  cnd = norm(double(inv(vpa(A))))*norm(A);
