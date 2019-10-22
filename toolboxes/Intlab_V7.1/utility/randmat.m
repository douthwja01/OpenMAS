function A = randmat(n,condition,Pell)
%RANDMAT      (ill-conditioned) random matrix, coefficients in [-1,1]
%
%   A = randmat(n,condition)
%
% condition, if specified, is approximate condition of generated matrix
% dimension n must be greater than 1
%
% Extremely ill-conditioned matrices with integer entries generated following
%   S.M. Rump: A Class of Arbitrarily Ill-conditioned Floating-Point Matrices, 
%     SIAM Journal on Matrix Analysis and Applications (SIMAX, 12(4):645-653, 1991.
% or using unit lower triangular L and A'*A (until condition 1e100
%
%The call
%   A = randmat(n,condition,Pell)
%with Pell=1 forces to use construction by Pell's equation.
%
%Examples of ill-conditioned matrices include:
%
%cond_inf(A4) ~ 1.4e65
%   A4 = [ -5046135670319638,  -3871391041510136, -5206336348183639,  -6745986988231149 ; ...
%           -640032173419322,   8694411469684959,  -564323984386760,  -2807912511823001 ; ... 
%         -16935782447203334, -18752427538303772, -8188807358110413, -14820968618548534 ; ...
%          -1069537498856711, -14079150289610606,  7074216604373039,   7257960283978710 ];
%      
%cond_inf(A6) ~ 1.5e94     
%   A6 = [ 1810371096161830, -2342429902418850,  2483279165876947, -3747279334964253, -3262424857701958, -3083142983440029 ; ...
%         -4670543938170405, -1397606024304777,    60011496034489,  1689416277492541, -1500903035774466,  3966198838365752 ; ...
%         -1064600276464066, -7561599142730362,  4805299117903646, -6663945806445749, -7071076225762059,   -52156788818356 ; ...
%         13002500911063530,  2223055422646289, -1553584081743862, -5252100317685933,  7524433713832350, -6396043249912008 ; ...
%          -395183142691090, -2180846347388541,  1450541682835654, -3629498209925700, -1866168768872108,  1230298410784196 ; ...
%          2337744233608461,  1359019382927754,  1241733688092475,  1803080888028433, -2648685047371017, -7046414836143443 ];
%

% written  04/27/97     S.M. Rump
% modified 04/11/99     S.M. Rump  exact match of condition
% modified 04/02/02     S.M. Rump  dimension greater than 1 checked
% modified 10/27/03     S.M. Rump  extremely ill-conditioned matrices added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 02/23/05     S.M. Rump  parameter Pell added
% modified 03/11/05     S.M. Rump  Samples A4, A6 added
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 03/01/08     S.M. Rump  condition numbers near overflow
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 01/28/11     S.M. Rump  threshhold for extremely ill-co now 1e16
% modified 05/28/11     S.M. Rump  rectangular matrices up cond 1e16
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==2
    Pell = 0;
  end
  A = 1e200;
  
  if length(n)>1
    m = n(1);
    n = n(2);
    if condition>1e16
      error('rectangular matrices only up to condition number 1e16')
    end
  else
    m = n;
  end

  if condition<=1e16
    
    % the usual generation of matrices with anticipated condition via svd
    if n==1
      error('dimension for randmat must be greater than 1')
    end
    
    if nargin==1
      A = random(n);
    else
      s = exp( [ 0 rand(1,min(m,n)-2) 1 ] * log(condition) ) / condition;
      if m>n
        A = randorth(m)*[diag(s);zeros(m-n,n)]*randorth(n);
      else
        A = randorth(m)*[diag(s) zeros(m,n-m)]*randorth(n);
      end
    end
    
  elseif ~Pell
    
    if condition <= 1e111            
      % generation via unit lower triangular L and A'*A, works until approximately condition 1e100
      % singular values well distributed between smallest and largest
      
      log10cond = log10(condition);
      k = 0;
      while log10cond > 14
        k = k+1;
        log10cond = log10cond/2;
      end
      log10c = 0;
      e = 1;
      I = eye(n);
      index = 0;
      while abs(log10cond-log10c)>1
        index = index+1;
        if index>5
          index = 0;
          if log10c<log10cond
            e = 0.5 + (e-0.5)*1.1;
          else
            e = 0.5 + (e-0.5)*0.9;
          end
        end
        A = tril(round(random(n,{-e,e})),-1) + I;
        L = tril(round(random(n,{-e,e})),-1) + I;
        A = L'*A;
        log10c = log10(cond(A));
      end
      for i=1:k
        A = A*A';
      end
%       A = L*A;
      N = n;

    end
  end

  if Pell | ( max(abs(A(:)))>1e16 )

    % generation via Pell's equation, only one very small singular value
    % generation of very ill-conditioned matrices with integer elements and determinant +/- 1
    % Careful: final dimension may be larger than wanted
    % final condition number in general larger than anticipated

    % M should not be too large, otherwise final matrix will not be exactly representable
    M = 1e14;                    % maximum size of individual component of PP, QQ

    while 1
      % generate random k for Pell's equation
      k = randint(19);
      succ = 0;
      if sqr(sqrt(k))~=k
        % generate initial solution for Pell's equation
        for Q=1:k^2
          P = round(sqrt(k*Q^2));
          if P^2-k*Q^2==1
            succ = 1;
            break
          end
        end
        if succ
          break
        end
      end
    end

    beta = randint(20)+1;       % generate basis for P,Q randomly
    P0 = normalize(P,beta);
    Q0 = normalize(Q,beta);
    kQ0 = k*Q0;

    % generate P,Q with P^2 of size condition
    %  for really large condition numbers, use
    %    while 2*length(P)*log(beta)<condition*log(10)
    % in which case  10^condition is the anticipated condition number
    while 2*length(P)*log(beta)<log(condition)
      N = add( mul(P,P0,beta) , mul(Q,kQ0,beta) , beta );
      Q = add( mul(Q,P0,beta) , mul(P,Q0,beta) , beta );
      P = N;
    end

    % adapt length of P and Q
    lenP = length(P);
    lenQ = length(Q);
    LPQ = max(lenP,lenQ);
    P = [P zeros(1,LPQ-lenP)];
    Q = [Q zeros(1,LPQ-lenQ)];

    m = floor(n/2);             % anticipated length of PP and QQ
    m = max( m , ceil(LPQ*log(beta)/log(M)) );  % actual length for maximum components M
    cont = 1;
    while cont                    % try this value of m
      d = ceil(LPQ/m);            % digits per component
      sigma = beta^d;             % basis for PP, QQ
      PP = zeros(1,m);
      QQ = zeros(1,m);
      for i=1:m
        v =  1+(i-1)*d : min(i*d,LPQ) ;
        PP(i) = sum( P(v) .* beta.^(0:length(v)-1) );
        QQ(i) = sum( Q(v) .* beta.^(0:length(v)-1) );
      end

      % calculate matrix
      PP = fliplr(PP);
      QQ = fliplr(QQ);
      N = max(n,2*m);             % actual dimension of matrix
      A = zeros(N);
      A(1,1:m) = PP;
      A(2,1:m) = QQ;
      A(1,m+1:2*m) = k*QQ;
      A(2,m+1:2*m) = PP;
      for i=1:m-1
        A(i+2,i) = 1;
        A(i+2,i+1) = -sigma;
        A(m+1+i,m+i) = 1;
        A(m+1+i,m+i+1) = -sigma;
      end
      if odd(N), A(N,N) = 1; end

      % extra row and column operations for 'random' structure
      wng = warning;
      warning off
      index = 0;
      while cont & ( index<10 )
        index = index+1;
        L = tril(randint(-1,N),-1) + eye(N);
        Linv = round(inv(L));
        % make sure Linv = L^-1
        setround(-1), LLinv1 = L*Linv;
        setround(1), LLinv2 = L*Linv;
        setround(0)
        if  ( ~isequal(LLinv1,LLinv2) ) | ( ~isequal(LLinv1,eye(N)) ), continue, end
        U = triu(randint(-1,N),1) + eye(N);
        Uinv = round(inv(U));
        % make sure Uinv = U^-1
        setround(-1), UUinv1 = U*Uinv;
        setround(1), UUinv2 = U*Uinv;
        setround(0)
        if ( ~isequal(UUinv1,UUinv2) )  | ( ~isequal(UUinv1,eye(N)) ), continue, end
        B = L*U*A*L'*U';
        % make sure no rounding error in the computation of B
        setround(-1), B = L*U*A*L'*U';
        setround(1), B2 = L*U*A*L'*U';
        setround(0)
        cont = ~isequal(B,B2);
      end
      if cont                     % transformed matrix not exactly representable
        m = m+1;                  % increase m and therefore size of matrix
      else
        A = B;
        warning(wng)
      end
    end
  end
  
  if rndold
    setround(rndold)
  end

    
  
  % A little multiple precision arithmetic
  function res = add(A,B,beta)
  % long addition to base beta for positive A,B
  lenA = length(A);
  lenB = length(B);
  res = normalize( [ A zeros(1,lenB-lenA) ] + [ B zeros(1,lenA-lenB) ] , beta );
  
  
  function res = mul(A,B,beta)
  % long multiplication to base beta for positive A,B
  res = normalize( conv(A,B) , beta );
  
  
  function A = normalize(A,beta)
  % normalitation to base beta for positive A
  while any( A>=beta )
    q = floor( A/beta);
    A = A - q*beta;
    A(2:end) = A(2:end) + q(1:end-1);
    if q(end)>0
      A = [ A q(end) ];
    end
  end
  
