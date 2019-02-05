function stdfctsdata
%STDFCTSDATA  Generation of data for rigorous intval standard functions
%
%   stdfctsdata
%
%If relative errors are too large, function will stop prematurely
%and corresponding flag is stored in stdfctsdata.mat
%

% written  12/30/98     S.M. Rump
% modified 09/25/99     S.M. Rump  generated data file appended with
%                                  Matlab version, tanh added and major
%                                  revision, length of arrays <=16384
%                                  to support Matlab student version
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 03/09/02     S.M. Rump  Windows for progress etc. deleted (Matlab problems)
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 08/21/04     S.M. Rump  file naming to cure Matlab/Unix bug
% modified 06/08/07     S.M. Rump  EPSMAX and NoSuccess (thanks to Carlos Lopez)
% modified 09/10/07     S.M. Rump  INTLAB path, file name
% modified 10/14/08     S.M. Rump  allow routine to be interrupted (thanks to 
%                                    Zhang Jun Jie for pointing to this)
% modified 05/29/09     S.M. Rump  omit display of INTLABPATH
% modified 08/26/12     S.M. Rump  global variables removed, see
% modified 10/03/12     S.M. Rump  finished comment
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
% modified 12/03/12     S.M. Rump  INTLABdata reorganization
%

  setround(0)                           % set rounding to nearest

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  INTLAB_STDFCTS_EPSMAX = getappdata(0,'INTLAB_STDFCTS_EPSMAX');
  
  ALLDATA = [];                         % struct array for all data
  
  INTLAB_INTLABPATH = getappdata(0,'INTLAB_INTLABPATH');
  ALLDATA = storedata(ALLDATA,INTLAB_INTLABPATH);

  INTLAB_ROUNDING_TRYMODE = getappdata(0,'INTLAB_ROUNDING_TRYMODE');
  ALLDATA = storedata(ALLDATA,INTLAB_ROUNDING_TRYMODE);

  INTLAB_ROUNDING_PARAM = getappdata(0,'INTLAB_ROUNDING_PARAM');
  ALLDATA = storedata(ALLDATA,INTLAB_ROUNDING_PARAM);

  INTLAB_ENV_VAR_BLAS_VERSION = getappdata(0,'INTLAB_ENV_VAR_BLAS_VERSION');
  ALLDATA = storedata(ALLDATA,INTLAB_ENV_VAR_BLAS_VERSION);

  INTLAB_INTVAL_POWER10 = getappdata(0,'INTLAB_INTVAL_POWER10');
  ALLDATA = storedata(ALLDATA,INTLAB_INTVAL_POWER10);
  

  % Maximum allowed value of relative error of built-in standard functions
  % for test set. If error is larger, stdfctsdata ends without success,
  % i.e. rigorous standard functions cannot be used on that computer
  setappdata(0,'INTLAB_STDFCTS_EPSMAX',1e-15);

  % save long environment
  errorterm = longinit('ErrorTerm',0);
  longinit('WithErrorTerm',0)
  logbetaold = INTLAB_LONG_LOGBETA;

  % initalize long toolbox
  INTLAB_LONG_LOGBETA = 25;
  INTLAB_LONG_BETA = 2^25;
  prec = longprecision;
  longprecision(35);


  disp('Finished about 0 %')
  
  disp('Check special function values ')
  success = CheckSpecialValues
  if ~success
    NoSuccess(ALLDATA,errorterm,prec,logbetaold);
  end
  disp(' ')
  disp('Finished about 2 %')
  
  disp('Error estimates to be computed:  Log ')
  [success,ALLDATA] = EstimateLog(ALLDATA);
  if ~success
    NoSuccess(ALLDATA,errorterm,prec,logbetaold);
  end
  disp(' ')

  disp('Error estimates to be computed:  Atan ')
  [success,ALLDATA] = EstimateAtan(ALLDATA);
  if ~success
    NoSuccess(ALLDATA,errorterm,prec,logbetaold);
  end
  disp(' ')
 
  disp('Error estimates to be computed:  ExpSinCosTanSinhTanh ')
  [success,ALLDATA] = EstimateExpSinCosTanSinhTanh(ALLDATA);
  if ~success
    NoSuccess(ALLDATA,errorterm,prec,logbetaold);
  end

  % error estimates successfully completed
  INTLAB_STDFCTS_SUCCESS = 1;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SUCCESS);
  
  % store data in file
  savedata(ALLDATA);
  disp('Error estimates in stdfctsdata successfully completed')
  
  

function NoSuccess(ALLDATA,errorterm,prec,logbetaold)
% quit stdfctsdata with no success

  INTLAB_STDFCTS_SUCCESS = 0;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SUCCESS);
  
  % make sure all constants not usable
  INTLAB_STDFCTS_LOG.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_LOG);
  INTLAB_STDFCTS_ATAN.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_ATAN);
  INTLAB_STDFCTS_EXP.POW = NaN;
  INTLAB_STDFCTS_EXP.POWINF = NaN;
  INTLAB_STDFCTS_EXP.POWSUP = NaN;
  INTLAB_STDFCTS_EXP.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_EXP);
  INTLAB_STDFCTS_SIN.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SIN);
  INTLAB_STDFCTS_COS.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_COS);
  INTLAB_STDFCTS_TAN.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_TAN);
  INTLAB_STDFCTS_SINH.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SINH);
  INTLAB_STDFCTS_TANH.EPS = NaN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_TANH);
  
  % restore long environment
  longinit(errorterm,0)
  longprecision(prec);
  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');
  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');
  INTLAB_LONG_LOGBETA = logbetaold;
  INTLAB_LONG_BETA = 2^logbetaold;
  
  % save NaN data
  savedata(ALLDATA);

  error('RELATIVE ERRORS TOO LARGE; USE OF RIGOROUS STANDARD FUNCTIONS NOT POSSIBLE')




function success = CheckSpecialValues
%Check special values for built-in standard functions:
%
%   exp(0) = 1
%   log(1) = 0
%   sin(0) = 0
%   cos(0) = 1
%   tan(0) = 0
%   atan(0) = 0
%   sinh(0) = 0
%   tanh(0) = 0
%

  success = 1;
  success = success & CheckSpecialValues_('exp',0,1);
  success = success & CheckSpecialValues_('log',1,0);
  success = success & CheckSpecialValues_('sin',0,0);
  success = success & CheckSpecialValues_('cos',0,1);
  success = success & CheckSpecialValues_('tan',0,0);
  success = success & CheckSpecialValues_('atan',0,0);
  success = success & CheckSpecialValues_('sinh',0,0);
  success = success & CheckSpecialValues_('tanh',0,0);



function success = CheckSpecialValues_(fct,x,y)
%Check special values for standard functions:  y = fct(x)
%

  setround(0)
  success = ( feval(fct,x) == y );



function [success,ALLDATA] = EstimateAtan(ALLDATA)
%For 15-bit mantissa numbers  x = ( 1 : 2^15-1)/2^13,
%the following is computed:
%
% ATAN.EPS    update for true bounds, i.e. for atanx:=atan(x) in
%             rounding to nearest,
%               atanx - ATAN.EPS*atanx <= atan(x) <= atanx + ATAN.EPS*atanx
%

  setround(0)
  ATAN.EPS = 0;

  for i=0:31          % range ( 1 .. 2^15-1 ) / 2^13
    if i==0
      v = pow2( ( 1 : 2^10-1 )' , -13 );
    else
      v = pow2( ( i*2^10 : (i+1)*2^10-1 )' , -13 );
    end
    [ ATAN , success ] = EstimateAtan_(ATAN,v);
    disp(['Finished about ' int2str(26+(i+1)*1.5)  ' %'])
    if ~success
      return
    end
  end

  INTLAB_STDFCTS_ATAN = ATAN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_ATAN);



function [ ATAN , success ] = EstimateAtan_(ATAN,v)
%Estimation for input arguments in v
%

  setround(0)

  ys = atan(v);
  zs = ys + ( v - longtan(ys) )./(1+v.*v);   % correction for atan
  e = mag(long2intval( longtan(zs) - v ));   % abs(zs-atan(v)) <= 3*abs(e)
  Y = addlongerror(zs,3*e);                  % atan(v) in Y

  [ ATAN , success ] = checkerror('atan',v,Y,ATAN);



function y = longtan(x)
% slow but rigorous tan for long x to current long precision
% (only for atan correction)

  prec = longprecision;
  T = long(x);
  S = T;
  C = 1;
  i = 1;
  while 1
    i = i+1;
    T = T*x/i;
    switch mod(i,4)
      case 0,  C = C + T;
      case 1,  S = S + T;
      case 2,  C = C - T;
      case 3,  S = S - T;
    end
    if all( mag(long2intval(T))<10^(-prec) )
      C = addlongerror( C , 10^(-prec) );
      S = addlongerror( S , 10^(-prec) );
      break
    end
  end
  y = S/C;




function [success,ALLDATA] = EstimateLog(ALLDATA)
%For 14-bit mantissa numbers  x = ( 2^13 : 2^14-1)*2^(E-14), 0<=E<=1,
%the following is computed:
%
% LOG.EPS     update for true bounds, i.e. for logx:=log(x) in
%             rounding to nearest,
%               logx - LOG.EPS*abs(logx) <= log(x) <= logx + LOG.EPS*abs(logx)
%

  setround(0)
  LOG.EPS = 0;

  for E=0:1
    for i=0:7
      if ( E==1 ) & ( i==0 )    % avoid v==1
        v = pow2( ( 2^13+1 : 2^13+2^10-1 )' , E-14 );
      else
        v = pow2( ( 2^13+i*2^10 : 2^13+(i+1)*2^10-1 )' , E-14 );
      end
      [ LOG , success ] = EstimateLog_(LOG,v);
      disp(['Finished about ' int2str(2+8*E+(i+1)) ' %']);
      if ~success
        return
      end
    end
  end

  INTLAB_STDFCTS_LOG = LOG;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_LOG);




function [ LOG , success ] = EstimateLog_(LOG,v)
%Estimation for input arguments in v
%

  setround(0)
  ys = log(v);
  zs = ys + ( v*exp(long(-ys)) - 1 );    % long correction for log
  e = long2intval( exp(zs) - v );
  e = mag( 2*e./(1-2*e) );               % abs(zs-log(v)) <= abs(e)
  Y = addlongerror(zs,e);                % log(v) in Y

  [ LOG , success ] = checkerror('log',v,Y,LOG);




function [success,ALLDATA] = EstimateExpSinCosTanSinhTanh(ALLDATA)
%For 14-bit mantissa numbers  x = +/- (1:2^14-1)*2^(-14) and
%   x = + (1:2^14-1)*2^(e-14)  where 1<=e<=3  the following is computed:
%
% EXP.EPS     update for true bounds, i.e. for expx:=exp(x) in
%             rounding to nearest it is
%               expxs*(1-EXP.EPS) <= exp(x) <= expxs*(1+EXP.EPS)
%
%Furthermore, for -744 <= index <= 709,
%
%    EXP.POW(745+index) + EXP.POWINF(745+index) <= exp(745+index) <=
%    EXP.POW(745+index) + EXP.POWSUP(745+index)
%
%
%For positive 14-bit mantissa numbers  x = (1:(2^14*pi/4))*2^(-14),
%which includes pi/4, the following is computed:
%
% SIN.EPS     update for true bounds, i.e. for sinx:=sin(x) in
%             rounding to nearest it is
%               sinx*(1-SIN.EPS) <= sin(x) <= sinx*(1+SIN.EPS)
%
% and similarly for COS and TAN
%
%
%For positive 14-bit mantissa numbers  x = (1:2^14-1)*2^(E-14), 0<=E<=3
%the following is computed:
%
% SINH.EPS     update for true bounds, i.e. for sinhx:=sinh(x) in
%              rounding to nearest it is
%               sinhx*(1-SINH.EPS) <= sinh(x) <= sinhx*(1+SINH.EPS)
%
%
%For positive 14-bit mantissa numbers  x = (1:2^14)*2^(-14)
%the following is computed:
%
% TANH.EPS     update for true bounds, i.e. for tanhx:=tanh(x) in
%              rounding to nearest it is
%               tanhx*(1-TANH.EPS) <= tanh(x) <= tanhx*(1+TANH.EPS)
%

  E_EXP = ExpPowers;

  E_EXP.EPS = 0;
  E_SIN.EPS = 0;
  E_COS.EPS = 0;
  E_TAN.EPS = 0;
  E_SINH.EPS = 0;
  E_TANH.EPS = 0;

  for i=0:15            % range ( 0 .. 2^14-1 ) / 2^14
    v = pow2( ( i*2^10 : (i+1)*2^10-1 )' , -14 );
    [E_EXP,E_SIN,E_COS,E_TAN,E_SINH,E_TANH,success] = ...
      EstimateExpSinCosTanSinhTanh_(E_EXP,E_SIN,E_COS,E_TAN,E_SINH,E_TANH,v);
    disp(['Finished about ' int2str(74+(i+1)*1.5) ' %'])
    if ~success
      return
    end
  end

  INTLAB_STDFCTS_EXP = E_EXP;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_EXP);

  INTLAB_STDFCTS_SIN = E_SIN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SIN);

  INTLAB_STDFCTS_COS = E_COS;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_COS);

  INTLAB_STDFCTS_TAN = E_TAN;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_TAN);

  INTLAB_STDFCTS_SINH = E_SINH;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_SINH);

  INTLAB_STDFCTS_TANH = E_TANH;
  ALLDATA = storedata(ALLDATA,INTLAB_STDFCTS_TANH);




function G_EXP = ExpPowers
% Compute EXP.POW
%

  precold = longprecision;
  prec = 50;
  longprecision(prec);
  setround(0)

  offset = 745;
  G_EXP.POW = zeros(1454,1);
  G_EXP.POWINF = zeros(1454,1);
  G_EXP.POWSUP = zeros(1454,1);

  % Bounds for powers of e
  G_EXP.POW(offset) = 1;
  G_EXP.POWINF(offset) = 0;
  G_EXP.POWSUP(offset) = 0;

  E = exp(long(1),prec);

  G_EXP.POW(offset+1) = long2dble(E);
  c = long2intval( E - G_EXP.POW(offset+1) );
  G_EXP.POWINF(offset+1) = c.inf;
  G_EXP.POWSUP(offset+1) = c.sup;

  EE = E;
  for index=2:709
    EE = EE*E;
    G_EXP.POW(offset+index) = long2dble(EE);
    c = long2intval( EE - G_EXP.POW(offset+index) );
    G_EXP.POWINF(offset+index) = c.inf;
    G_EXP.POWSUP(offset+index) = c.sup;
  end

  E = exp(long(-1),prec);

  G_EXP.POW(offset-1) = long2dble(E);
  c = long2intval( E - G_EXP.POW(offset-1) );
  G_EXP.POWINF(offset-1) = c.inf;
  G_EXP.POWSUP(offset-1) = c.sup;

  EE = E;
  for index=-2:-1:-744
    EE = EE*E;
    G_EXP.POW(offset+index) = long2dble(EE);
    c = long2intval( EE - G_EXP.POW(offset+index) );
    G_EXP.POWINF(offset+index) = c.inf;
    G_EXP.POWSUP(offset+index) = c.sup;
  end

  prec = precold;
  longprecision(prec);
  setround(0)



function [E_EXP,E_SIN,E_COS,E_TAN,E_SINH,E_TANH,success] = ...
  EstimateExpSinCosTanSinhTanh_(E_EXP,E_SIN,E_COS,E_TAN,E_SINH,E_TANH,v)

  % simultaneous computation of exp, sin, cos, tan, sinh, cosh, tanh
  prec = longprecision;
  X = long(v);
  LongExp = long(1) + X;
  LongExp_ = long(1) - X;
  LongSin = X;
  LongCos = long(1);
  LongSinh = X;
  LongCosh = long(1);
  Term = X;
  i = 1;
  while 1
    i = i+1;
    Term = Term * X / i;                           % Term = X^i/i!
    if all( mag(long2intval(Term))<10^(-prec) )
      LongExp = addlongerror( LongExp , 10^(-prec) );
      LongExp_ = addlongerror( LongExp_ , 10^(-prec) );
      LongSin = addlongerror( LongSin , 10^(-prec) );
      LongCos = addlongerror( LongCos , 10^(-prec) );
      LongSinh = addlongerror( LongSinh , 10^(-prec) );
      LongCosh = addlongerror( LongCosh , 10^(-prec) );
      break
    end
    LongExp = LongExp + Term;
    switch mod(i,4)
      case 0, LongCos = LongCos + Term;
      case 1, LongSin = LongSin + Term;
      case 2, LongCos = LongCos - Term;
      case 3, LongSin = LongSin - Term;
    end
    switch mod(i,2)
      case 0, LongCosh = LongCosh + Term;
              LongExp_ = LongExp_ + Term;
      case 1, LongSinh = LongSinh + Term;
              LongExp_ = LongExp_ - Term;
    end
  end

  if v(1)==0
    % special value f(0)
    LongExp(1) = long(1);
    LongExp_(1) = long(1);
    LongSin(1) = long(0);
    LongCos(1) = long(1);
    LongSinh(1) = long(0);
    LongCosh(1) = long(1);
  end

  % exp(v) in LongExp, exp(-v) in LongExp_
  % sin(v) in LongSin, cos(v) in LongCos
  % sinh(v) in LongSinh, cosh(v) in LongCosh

  % error bounds for exp, positive arguments 0<=x<1
  [ E_EXP , success ] = checkerror('exp',v,LongExp,E_EXP);
  if ~success
    return
  end

  % error bounds for exp, negative arguments -1<x<=0
  [ E_EXP , success ] = checkerror('exp',-v,LongExp_,E_EXP);
  if ~success
    return
  end

  % error bounds for exp, positive arguments 1<=x<2
  LongExp1 = LongExp.*LongExp;
  [ E_EXP , success ] = checkerror('exp',2*v,LongExp1,E_EXP);
  if ~success
    return
  end

  % error bounds for exp, positive arguments 2<=x<4
  LongExp1 = LongExp1.*LongExp1;
  [ E_EXP , success ] = checkerror('exp',4*v,LongExp1,E_EXP);
  if ~success
    return
  end

  % error bounds for exp, positive arguments 4<=x<8
  LongExp1 = LongExp1.*LongExp1;
  [ E_EXP , success ] = checkerror('exp',8*v,LongExp1,E_EXP);
  if ~success
    return
  end

  index = ( v>.785398164 );     % exclude components >pi/4
  if any(index)                 % v>pi/4 nonempty

    w = v(~index);              % w = v(v<=pi/4)
    if ~isempty(w)              % isempty(w) <=> all(v>pi/4)

      lw = length(w);

      % error bounds for sin
      LongSin = LongSin(1:lw);
      [ E_SIN , success ] = checkerror('sin',w,LongSin,E_SIN);
      if ~success
        return
      end

      % error bounds for cos
      LongCos = LongCos(1:lw);
      [ E_COS , success ] = checkerror('cos',w,LongCos,E_COS);
      if ~success
        return
      end

      LongTan = LongSin/LongCos;
      [ E_TAN , success ] = checkerror('tan',w,LongTan,E_TAN);
      if ~success
        return
      end

    end

  else                          % all(v<pi/4)

    % error bounds for sin
    [ E_SIN , success ] = checkerror('sin',v,LongSin,E_SIN);
    if ~success
      return
    end

    % error bounds for cos
    [ E_COS , success ] = checkerror('cos',v,LongCos,E_COS);
    if ~success
      return
    end

    LongTan = LongSin/LongCos;
    [ E_TAN , success ] = checkerror('tan',v,LongTan,E_TAN);
    if ~success
      return
    end

  end

  % error bounds for tanh
  LongTanh = LongSinh/LongCosh;
  [ E_TANH , success ] = checkerror('tanh',v,LongTanh,E_TANH);
  if ~success
    return
  end

  % error bounds for sinh, 0<=x<1
  [ E_SINH , success ] = checkerror('sinh',v,LongSinh,E_SINH);
  if ~success
    return
  end

  % error bounds for sinh, 1<=x<2
  Sh2 = 2*LongSinh*LongCosh;
  LongCosh = LongSinh*LongSinh + LongCosh*LongCosh;
  LongSinh = Sh2;
  [ E_SINH , success ] = checkerror('sinh',2*v,LongSinh,E_SINH);
  if ~success
    return
  end

  % error bounds for sinh, 2<=x<4
  Sh2 = 2*LongSinh*LongCosh;
  LongCosh = LongSinh*LongSinh + LongCosh*LongCosh;
  LongSinh = Sh2;
  [ E_SINH , success ] = checkerror('sinh',4*v,LongSinh,E_SINH);
  if ~success
    return
  end

  % error bounds for sinh, 4<=x<8
  LongSinh = 2*LongSinh*LongCosh;
  [ E_SINH , success ] = checkerror('sinh',8*v,LongSinh,E_SINH);
  if ~success
    return
  end




function [ INTLAB_DATA , success ] = checkerror(f,v,Ylong,INTLAB_DATA)
% check floating point error against long, check max error

  INTLAB_STDFCTS_EPSMAX = getappdata(0,'INTLAB_STDFCTS_EPSMAX');
  success = 1;

  % fl-pt approximation
  setround(0)
  cflpt = feval(f,v);

  D = long2intval( Ylong-cflpt );        % error cflpt against f(x)
  index = ( D.inf<0 );                   % should cflpt<=Ylong
  if any(index)
    setround(1)
    E = mag(D(index)) ./ abs(cflpt(index));
    % compute error for lower bound
    INTLAB_DATA.EPS = max( INTLAB_DATA.EPS , max(E) );
    if INTLAB_DATA.EPS > INTLAB_STDFCTS_EPSMAX
      success = 0;
      return
    end
  end

  D = long2intval( cflpt-Ylong );        % error cflpt against f(x)
  index = ( D.inf<0 );                   % should Ylong<=cflpt
  if any(index)
    setround(1)
    E = mag(D(index)) ./ abs(cflpt(index));
    % compute error for upper bound
    INTLAB_DATA.EPS = max( INTLAB_DATA.EPS , max(E) );
    if INTLAB_DATA.EPS > INTLAB_STDFCTS_EPSMAX
      success = 0;
      return
    end
  end

  setround(0)



function ALLDATA = storedata(ALLDATA,DATA,see)
% saves dataname and value into ALLDATA

  datanew.name = inputname(2);
  datanew.value = DATA;
  ALLDATA = [ALLDATA datanew];
  if nargin<3
    see = 0;
  end
  if see
    disp(datanew.name)
    datanew.value
  end

  
function savedata(ALLDATA)
% write data into file

  fname = intvalinit('INTLABdata');
  evalstr = [ 'save ' '''' fname '''' ];   % make sure to allow blanks in fname
  for i=1:length(ALLDATA)
    if ischar(ALLDATA(i).value)      % take care of INTLABPATH
      eval([ ALLDATA(i).name ' = ''' ALLDATA(i).value ''';'])
    else
      eval([ ALLDATA(i).name ' = ALLDATA(i).value ;'])
    end
    evalstr = [ evalstr ' ' ALLDATA(i).name ];
  end
  eval(evalstr);
  disp(' ')
  disp([ '===> Data written to ' fname ])
  disp(' ')
