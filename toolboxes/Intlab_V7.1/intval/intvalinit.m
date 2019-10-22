function res = intvalinit(param,see,INTLABPATH)
%INTVALINIT   Initialization and defaults for intval toolbox
%
%   intvalinit(param,see)
%
%possible values for param:
%
%  'Reference'               Proper reference to INTLAB
%
%  'License'                 Gives license information
%
%  'Version'                 Gives INTLAB version information
%
%  'Init'                    Initialize INTLAB constants (for startup file)
%
%  'CheckRounding'           Check directed rounding
%
%  'StoreDefaultSetting'     Store current setting as default setting, see
%                              also "intlabsetting"
%
%Default display of intervals:
%  'Display_'                Display with uncertainty (e.g. 3.14_) but change to inf/sup for real
%                              or to mid/rad for complex input if intervals too wide 
%  'Display__'               Display with uncertainty (e.g. 3.14_) independent of width of input
%  'DisplayInfsup'           Display infimum/supremum (e.g. [ 3.14 , 3.15 ])
%  'DisplayMidrad'           Display midpoint/radius (e.g. < 3.14 , 0.01 >)
%  'Display'                 res = 'Display_'
%                                  'Display__'
%                                  'DisplayInfsup'
%                                  'DisplayMidrad'
%
%Default exception handling for real interval standard functions:
%  'RealStdFctsExcptnAuto'   Complex interval stdfct used automatically for
%                                real interval input out of range (without warning)
%  'RealStdFctsExcptnWarn'   As 'RealStdFctsExcptnAuto', but with warning
%  'RealStdFctsExcptnNaN'    Result NaN for real interval input out of range
%  'RealStdFctsExcptnIgnore' Input arguments out of range ignored
%  'RealStdFctsExcptn'       res = 'RealStdFctsExcptnAuto'
%                                  'RealStdFctsExcptnWarn'
%                                  'RealStdFctsExcptnNaN'
%                                  'RealStdFctsExcptnIgnore'
%  'RealStdFctsExcptnOccurred'   Flag: result 1 iff input argument was out of range 
%                                while option 'RealStdFctsExcptnIgnore' was set; 
%                                Flag is reset after check
%
%Default linear system solver with thin input matrix:
%  'DoubleResidual'          Residual in double precision
%  'ImprovedResidual'        Simulated higher precision residual improvement for improved
%                              accuracy, for details see intval\lssresidual
%  'QuadrupleResidual'       Quadruple precision residual, see accsumdot\Dot_
%  'Residual'                res = 'DoubleResidual'
%                                  'ImprovedResidual'
%                                  'QuadrupleResidual'
%
%Default thick real interval times thick real interval:
%  'FastIVmult'              Fast algorithm by midpoint/radius arithmetic
%                              absolutely rigorous, worst case overestimation 1.5
%  'SharpIVmult'             Slow algorithm by infimum/supremum arithmetic
%                              slow for larger matrices due to interpretation overhead
%  'IVmult'                  res = 'FastIVmult'
%                                  'SharpIVmult'
%
%Default warning handling for x(i)=c, where x real interval array and c complex
%  'RealComplexAssignAuto'   Real interval arrays automatically transformed to
%                                complex array if a component is assigned a
%                                complex value (without warning)
%  'RealComplexAssignWarn'   As 'RealComplexAssignAuto', but with warning
%  'RealComplexAssignError'  Assignment of a complex value to a component of a
%                                real interval array causes an error
%  'RealComplexAssign'       res = 'RealComplexAssignAuto'
%                                  'RealComplexAssignWarn'
%                                  'RealComplexAssignError'
%Default warning handling for complex interval defined by infsup(zinf,zsup)
%  'ComplexInfSupAssignWarn'     Warning issued if complex interval defined by infsup(zinf,zsup)
%  'ComplexInfSupAssignNoWarn'   No warning issued if complex interval defined by infsup(zinf,zsup)
%  'ComplexInfSupAssign'         res = 'ComplexInfSupAssignWarn'
%                                      'ComplexInfSupAssignWarn'
%
%A corresponding message is printed when changing a default mode unless
%  it is explicitly suppressed with the (optional) parameter see=0.
%

% written  10/16/98     S.M. Rump
% modified 10/24/99     S.M. Rump  adapted to INTLAB V3, data files
%                                  appended with Matlab version, welcome
% modified 09/02/00     S.M. Rump  real/complex asignment
%                                  rounding unchanged after use
%                                  std fcts switch non-rigorous corrected
%                                    (thanks to S. Christiansen and N. Albertsen)
% modified 11/16/01     S.M. Rump  isieee removed
% modified 12/15/01     S.M. Rump  displaywidth added
% modified 02/15/02     S.M. Rump  Execution terminated if rounding fails
%                                  Automatic choice of setround .m / .dll
%                                  system_dependent('setround',0.5) for round to nearest (thanks to J.A. van de Griend)
% modified 03/09/02     S.M. Rump  Windows for stdfctsdata etc. deleted (Matlab problems)
% modified 04/02/02     S.M. Rump  Residual switch added
% modified 10/09/02     S.M. Rump  Rounding switch by global variable
% modified 11/30/03     S.M. Rump  second stage and quadruple prec. residual added
% modified 01/10/04     S.M. Rump  default display changed to infsup
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    InfNanFlag added
% modified 08/15/04     S.M. Rump  various little changes
% modified 08/21/04     S.M. Rump  INTLAB version added, file naming to cure Matlab/Unix problems
%                                    (thanks to George Corliss and Annie Yuk for pointing to this)
% modified 10/04/04     S.M. Rump  Warning for infsup(zinf,zsup) added
% modified 01/06/05     S.M. Rump  Check for Atlas BLAS added
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/02/05     S.M. Rump  'realstdfctsexcptnignore' added
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed wording (thanks to Jiri Rohn)
% modified 06/30/06     S.M. Rump  Sequence of checkrounding changed
% modified 12/03/06     S.M. Rump  Sparse Bug global flag (thanks to Arnold)
% modified 12/05/06     S.M. Rump  flag Display__, helpp file added
% modified 02/18/07     S.M. Rump  welcome, thanks to J. Rohn and C. Linsenmann
% modified 09/06/07     S.M. Rump  approximate std fcts removed,
%                                    version number, 
%                                    data files forced into INTLAB directory,
%					               	 help text for atlas routines changed
% modified 02/16/08     S.M. Rump  option 'references' added
%                                    .datenum in version detection (fixes R2008a bug)
%                                    path for rounding
% modified 05/06/08     S.M. Rump  comments to rounding adapted
%                                    warning off in CheckRounding
%                                    datenum removed
% modified 05/19/08     S.M. Rump  Reorganization of setround
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 09/28/08     S.M. Rump  check for INTLAB_SPARSE_BUG improved
% modified 10/14/08     S.M. Rump  INTLAB path included in stdfctsdatafile
% modified 10/19/08     S.M. Rump  NaN- and ignore-mode defined, RealStdFctsExcptn adapted
% modified 05/02/09     S.M. Rump  typo
% modified 05/29/09     S.M. Rump  'Reference' added
% modified 11/02/09     S.M. Rump  check for multi-core problems
% modified 11/21/09     S.M. Rump  check for R2009b
% modified 02/23/09     S.M. Rump  reorganization of rounding check (thanks to Viktor Härter)
%                                  references reorganized
% modified 04/07/10     S.M. Rump  check for BLAS_VERSION to avoid problems with omp
% modified 06/18/10     S.M. Rump  wording in comments
% modified 08/07/10     S.M. Rump  upper case Dot_
% modified 12/06/10     S.M. Rump  feature_omp removed on 64-bit Unix systems
%                                     (thanks to Frank Schmidt, Chemnitz)
% modified 02/17/11     S.M. Rump  pause replaced by input('...'), cures Matlab bug
% modified 02/20/11     S.M. Rump  switching rounding in multi-core (Thanks to Dr. Ozaki)
% modified 04/03/12     S.M. Rump  typo
% modified 08/24/12     S.M. Rump  SECOND_STAGE superflous due to redesign of verifylss
% modified 08/26/12     S.M. Rump  INTLABPATH corrected
% modified 08/26/12     S.M. Rump  global variables removed
% modified 09/07/12     S.M. Rump  initialize long package
% modified 10/03/12     S.M. Rump  intlabsetting and start pages
% modified 10/04/12     S.M. Rump  warning unaltered
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
% modified 10/16/12     S.M. Rump  Warning std fcts out of range
% modified 11/06/12     S.M. Rump  welcome text
% modified 11/12/12     S.M. Rump  check for fatal error
% modified 11/13/12     S.M. Rump  check for versions with flaws/bugs
% modified 11/14/12     S.M. Rump  extra check switching rounding back and forth
% modified 12/03/12     S.M. Rump  Reorganization of first start
% modified 12/07/12     S.M. Rump  Cure Unix behaviour in detecting correct rounding
% modified 12/09/12     S.M. Rump  Best rounding stored permanently, improved check
%                                    for Unix systems (thanks to M. Lange),
%                                    improved messages and very first start
% modified 02/28/12     S.M. Rump  Some comments corrected
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ( nargin==1 )
    see = 1;
  end

  switch lower(param)
      
  %%%%%%%%%%%%%%%%% test for directed rounding  %%%%%%%%%%%%%%%%%%
  case 'checkrounding'
    try
      wngNameConflict = warning('off','MATLAB:dispatcher:nameConflict');
      wngMexExt = warning('off','MATLAB:dispatcher:ShadowedMEXExtension');
    end
    %%% remove all setround directories from path
    p = [ pathsep path pathsep ];       % complete path
    k = findstr('setround',p);
    while ~isempty(k)
      k = findstr(pathsep,p(1:k(1)));
      p = p(k(end)+1:end);
      k = findstr(pathsep,p);
      p = p(1:k(1)-1);
      rmpath(p);
      p = [ pathsep path pathsep ];     % complete path
      k = findstr('setround',p);
    end

    % find setround path
    INTLAB_INTLABPATH = getappdata(0,'INTLAB_INTLABPATH');
    setround_path = [ INTLAB_INTLABPATH 'setround' filesep ];
   
    % Check for very first call
    fname = intvalinit('INTLABdata');
    VeryFirstCall = ~exist(fname,'file');
    
    if VeryFirstCall
      INTLAB_ROUNDING_TRYMODE = 0;
      intvalmessages('FindBestRounding')
      dummy = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% check rounding
      % setround_mex          fpu-flag by assembler
      % setround_sse          fpu- and sse-flag by assembler (files setround)
      % setround_mex_sse      fpu- and sse-flag by assembler (files setround_mex)
      % setround_mex_sse4     fpu- and sse-flag by assembler (files setround_mex4)
      % feature               fpu-flag by feature
      % feature_sse           fpu- and sse-flag by feature
      % feature_sse4          fpu- and sse-flag by feature multi-core
      % system_dependent      fpu-flag by system_dependent
      % system_dependent_sse  fpu- and sse-flag by system_dependent
      % system_dependent_sse4 fpu- and sse-flag by system_dependent multi-core
      % setround_omp          fpu- and sse-flag by fesetround (C99), multi-core
      % feature_omp           fpu-flag by feature and sse-flag as above, multi-core
      %
      % parameter 0.5 or 'nearest' for rounding to nearest by system_dependent
      INTLAB_ROUNDING_PARAM = { {'system_dependent',0.5};
        {'system_dependent','nearest'};
        {'setround_mex',dummy};
        {'setround_mex_sse',dummy};
        {'setround_sse',dummy};
        {'feature',dummy};
        {'feature_sse4','setround_mex_sse4',dummy};
        {'system_dependent_sse','setround_mex_sse',0.5};
        {'system_dependent_sse4','setround_mex_sse4',0.5};
        {'system_dependent_sse','setround_mex_sse','nearest'};
        {'system_dependent_sse4','setround_mex_sse4','nearest'}
        {'feature_sse','setround_mex_sse',dummy}; % might cause problems
        {'setround_omp',dummy};   % might cause problems
        {'feature_omp',dummy};    % do not use on 64-bit Unix system
        };      
      if isequal(computer,'GLNXA64')
        %       disp('===> On some 64-bit Unix systems a strange behaviour was observed:          ***')
        %       disp('===> The system "freezes" completely (thanks to Frank Schmidt in Chemnitz). ***')
        %       disp('===> If this happens, please comment line 226 and uncomment line 227        ***')
        %       disp('===> in intvalinit.m before restarting INTLAB.                              ***')
        %       input('Press enter to continue.')
        INTLAB_ROUNDING_PARAM(end) = [];      % remove feature_omp, causes problems under 64-bit Unix
      end
      if isempty(getenv('BLAS_VERSION'))
        INTLAB_ENV_VAR_BLAS_VERSION = 0;
      else
        INTLAB_ENV_VAR_BLAS_VERSION = 1;
        disp('%%%%%%%%%%%%% Installation of INTLAB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
        disp('%%% ')
        disp('%%% It seems the system variable "BLAS_VERSION" is still set.')
        disp('%%% This was necessary in older Matlab versions; now the ')
        disp('%%% Intel Math Kernel library is working with directed rounding')
        disp('%%% and is preferable because it supports multi-core and is much')
        disp('%%% faster.')
        disp('%%% Please unset the system variable "BLAS_VERSION" in your')
        disp('%%% operating system (sorry, can''t do that in INTLAB).')
        disp('%%%')
        disp('%%% Press Enter to continue.')
        input('%%%')
      end
      % very first initialization
      INTLAB_STDFCTS_SUCCESS = 0;
      save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
        'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
    else
      load(fname);
      % check whether previous try failed
      if ( INTLAB_ROUNDING_TRYMODE>0 ) & ( ~isinf(INTLAB_ROUNDING_TRYMODE) )
        % INTLAB_ROUNDING_TRYMODE    0   just starting, cannot happen here
        %                          ~=0   routine INTLAB_ROUNDING_TRYMODE failed
        %                          inf   selection process finished
        intvalmessages('BestRoundingTestFailed',INTLAB_ROUNDING_PARAM{INTLAB_ROUNDING_TRYMODE})
        % delete that rounding from INTLAB_ROUNDING_PARAM
        index = [ 1:INTLAB_ROUNDING_TRYMODE-1 INTLAB_ROUNDING_TRYMODE+1:length(INTLAB_ROUNDING_PARAM)];
        INTLAB_ROUNDING_PARAM = INTLAB_ROUNDING_PARAM(index);
        % intermediate state
        save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
          'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
      end
    end
  
    IMKL = ones(1,length(INTLAB_ROUNDING_PARAM)); % only IMKL failed
    T = IMKL;                           % timing for rounding procedures                                        
    mthr = IMKL;                        % rounding control for multi-core 
                                 
    if ~isinf(INTLAB_ROUNDING_TRYMODE)  % find best switching rounding routine
      if INTLAB_ROUNDING_TRYMODE==-1;
        disp('===> A previous call of INTLAB produced the following message')
        disp(' ')
        intvalmessages('MULTICOREerror');
        input(' ')
        exit;
      end
      if INTLAB_ROUNDING_TRYMODE==-2;
        disp('===> A previous call of INTLAB produced the following message')
        disp(' ')
        intvalmessages('IMKLerror');
        input(' ')
        exit;
      end
      if INTLAB_ROUNDING_TRYMODE==-3;
        disp('===> A previous call of INTLAB produced the following message')
        disp(' ')
        intvalmessages('ROUNDerror');
        input(' ')
        exit;
      end
      if INTLAB_ROUNDING_TRYMODE==-4;
        disp('===> A previous call of INTLAB produced the following message')
        disp(' ')
        intvalmessages('FatalError');
        input(' ')
        exit
      end
      for i=1:length(INTLAB_ROUNDING_PARAM)
        INTLAB_ROUNDING_TRYMODE = i;    % try INTLAB_ROUNDING_PARAM{i}
        % Save setting because ...
        save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
          'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
        % ... the next call might cause a hard stop
        [mthr(i),IMKL(i),T(i)] = trysetround(setround_path,INTLAB_ROUNDING_PARAM{i},i);
      end
      % Test succeeded, store environment
      INTLAB_ROUNDING_TRYMODE = inf;      % all roundings succeeded
      if any(mthr)              % 18.02.2010 VH
        mthr = logical(mthr); % mthr(i)==1 -> i-th round function is o.k.
        T(~mthr) = inf;       % in multi-core mode. If any rounding
        IMKL(~mthr) = inf;    % function correct at multi-core stage, take this fcn.
      else
        try
          wngThread = warning('off','MATLAB:maxNumCompThreads:Deprecated');
          maxNumCompThreads(1);   % else try to set the number of threads to one.
          warning(wngThread.state,'MATLAB:maxNumCompThreads:Deprecated');
        catch
          INTLAB_ROUNDING_TRYMODE = -1;
          save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
                     'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
          intvalmessages('MULTICOREerror');
          input(' ')
          exit;
        end
      end
      [minT,imin] = min(T);		          % choose fastest version
      failed = isequal(minT,inf);	      % no rounding setting working
      IMKL = max(IMKL);			          % IMKL=1: only IMKL failed
      if failed                           % rounding does not work properly
        if IMKL                           % error in rounding only using IMKL
          INTLAB_ROUNDING_TRYMODE = -2;
          save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
                     'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
          intvalmessages('IMKLerror');
          input(' ')
          exit;
        else
          INTLAB_ROUNDING_TRYMODE = -3;
          save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
                     'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
          intvalmessages('ROUNDerror');
          input(' ')
          exit
        end
      end
      INTLAB_ROUNDING_PARAM = INTLAB_ROUNDING_PARAM{imin};
      save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
        'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
    end
          
    % Fastest rounding found and stored; check again for safety
    for i=1:length(INTLAB_ROUNDING_PARAM)-1	    % restore fastest setting
      path( [ setround_path INTLAB_ROUNDING_PARAM{i} ] , path );
    end
    path(path);			    	    % make sure paths are correct
    INTLAB_ROUND_TO_NEAREST = INTLAB_ROUNDING_PARAM{end};
    setappdata(0,'INTLAB_ROUND_TO_NEAREST',INTLAB_ROUND_TO_NEAREST);
    failed = testrounding(0);        % last check for safety
    if failed                         % This should never happen
      INTLAB_ROUNDING_TRYMODE = -4;
      save(fname,'INTLAB_ROUNDING_PARAM','INTLAB_ROUNDING_TRYMODE',...
        'INTLAB_ENV_VAR_BLAS_VERSION','INTLAB_STDFCTS_SUCCESS')
      intvalmessages('FatalError');
      input(' ')
      exit
    end
      
    disp(' ')
    disp('===> rounding checked and no errors detected')
    warning(wngNameConflict.state,'MATLAB:dispatcher:nameConflict');
    warning(wngMexExt.state,'MATLAB:dispatcher:ShadowedMEXExtension');

  %%%%%%%%%%%%%%%%%%% license text message %%%%%%%%%%%%%%%%%%%%%
  case 'license'
      intvalmessages('LIC') 

  %%%%%%%%%%%%%%%%%% reference text message %%%%%%%%%%%%%%%%%%%%  
  case 'reference'
    intvalmessages('reference') 

  %%%%%%%%%%%%%%%%%%% store default setting %%%%%%%%%%%%%%%%%%%%
  case 'storedefaultsetting'
    setappdata(0,'INTLAB_SETTING',intlabsetting(0));
  
  %%%%%%%%%%%%%%%%%% set/get INTLAB version %%%%%%%%%%%%%%%%%%%%
  case 'version'
    INTLAB_VERSION = getappdata(0,'INTLAB_VERSION');
    if isempty(INTLAB_VERSION)
      INTLAB_INTLABPATH = getappdata(0,'INTLAB_INTLABPATH');
      in = fopen([INTLAB_INTLABPATH 'Contents.m'],'r'); 
      strContents = fscanf(in,'%c');    % the file Contents.m
      index = findstr('INTLAB_Version_',strContents); 
      index = index(end);               % the last mentioned Version
      str = strContents(index:index+25);
      i = findstr(' ',str);
      INTLAB_VERSION = str(16:i(1)-1);  % the version number
      setappdata(0,'INTLAB_VERSION',INTLAB_VERSION);
    end
    if see==0                           % internal call: only version number
      res = INTLAB_VERSION;
    else
      res = [ 'Intlab_Version_' INTLAB_VERSION ];
    end
    
  %%%%%%%%%%%%%%%%%%%%%% INTLAB variables %%%%%%%%%%%%%%%%%%%%%%%%
  case 'init'
      
    % Internal call with three parameters
    if nargin==3
      if ~isequal(INTLABPATH(end),filesep)
        INTLABPATH = [ INTLABPATH filesep ];
      end
      INTLAB_INTLABPATH = INTLABPATH;
      setappdata(0,'INTLAB_INTLABPATH',INTLAB_INTLABPATH);
    else
      INTLAB_INTLABPATH = getappdata(0,'INTLAB_INTLABPATH');
    end
    if see
      disp('===> Initialization of INTLAB toolbox')
    end

    %%% Check for very first call %%%
    VeryFirstCall = ~exist(intvalinit('INTLABdata'),'file');
    
    if VeryFirstCall
      versioncheck('');                 % check known bugs
      versioncheck('R2008b');           % rounding probably not working
      versioncheck('R2009a');           % rounding probably not working
    else
      % Always check for erroneous versions
      versioncheck('');                 % for safety: check various known bugs
    end
    
    intvalinit('checkrounding');        % check rounding
    setround(0)                         % initialize rounding mode
       
    load(intvalinit('INTLABdata'));
%     setappdata(0,'INTLAB_INTLABPATH',INTLAB_INTLABPATH);
    setappdata(0,'INTLAB_ROUNDING_TRYMODE',INTLAB_ROUNDING_TRYMODE);
    setappdata(0,'INTLAB_ROUNDING_PARAM',INTLAB_ROUNDING_PARAM);
    setappdata(0,'INTLAB_ENV_VAR_BLAS_VERSION',INTLAB_ENV_VAR_BLAS_VERSION);
    firstcall = ( ~INTLAB_STDFCTS_SUCCESS );
    
    setappdata(0,'INTLAB_INTVAL_DISPLAY','DisplayInfsup');
    
    INTLAB_SPARSE_BUG = ~any(any(isnan(full(speye(3)*inf))));
    setappdata(0,'INTLAB_SPARSE_BUG',INTLAB_SPARSE_BUG);
    
    % default width of display
    setappdata(0,'INTLAB_DISPLAY_WIDTH',120);
    
    % switch stdfct input out of real range
    setappdata(0,'INTLAB_STDFCTS_EXCPTN',1); % switch to complex interval stdfct with warning

    % flag: input out of range occurred while
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',0);   % 'RealStdFctsExcptnIgnore' active

    % computation of residuals
    setappdata(0,'INTLAB_INTVAL_RESIDUAL',1); % residual initialized to be more accurate
    
    % computation of interval matrix products
    setappdata(0,'INTLAB_INTVAL_IVMULT',0);   % IVmult initialized to be fast

    % Real interval array(i) = complex
    setappdata(0,'INTLAB_STDFCTS_RCASSIGN',1); % type conversion to complex (with warning)

    % largest x with exp(x) finite
    setappdata(0,'INTLAB_STDFCTS_LOGREALMAX',(12486629536300000+30718)*2^-44);

    %%% initialize constants for I/O %%%
    if firstcall
      disp(' ')
      disp('**** INTLAB: Generation of I/O data')
      INTLAB_INTVAL_POWER10 = initpower10;
      setappdata(0,'INTLAB_INTVAL_POWER10',INTLAB_INTVAL_POWER10);
      disp('**** INTLAB: I/O data successfully generated')
      %%% initialize constants for standard functions %%%
      intvalmessages('STDFCNmissed')
      longinit('init'); 	% make sure long package is initialized
      stdfctsdata
    end
   
    stdfctsinit
    INTLAB_STDFCTS_SUCCESS = getappdata(0,'INTLAB_STDFCTS_SUCCESS');
    if INTLAB_STDFCTS_SUCCESS
      if see
        disp('===> INTLAB constants for rigorous standard functions successfully computed')
      end
    else
      intvalmessages('CONSTerror')
      input(' ')
      exit
    end
    
    disp(' ')
    if firstcall
      intvalinit('Reference')
      disp('To see this reference, please use intvalinit(''Reference'').')
      disp(' ')
      input('Please press Enter to start INTLAB.')
      disp(' ')
      disp(' ')
    end
      
    disp(' ')
    v = intvalinit('version',0);
    v = [v blanks(7-length(v))];
    intvalmessages('WELCOME',v) 
    
  %%%%%%%%%%%%%% print display settings of INTVAL %%%%%%%%%%%%%%
  case 'display'
    res = getappdata(0,'INTLAB_INTVAL_DISPLAY');

  %%%%%%%%%%%%%%% set display settings for INTVAL %%%%%%%%%%%%%%
  case 'display_'
    setappdata(0,'INTLAB_INTVAL_DISPLAY','Display_');
    if see
      disp(sprintf(['===> Default display of intervals with uncertainty (e.g. 3.14_), changed \n' ...
                    '        to inf/sup or mid/rad if input too wide ']))
    end

  case 'display__'
    setappdata(0,'INTLAB_INTVAL_DISPLAY','Display__');
    if see
      disp(sprintf(['===> Default display of intervals with uncertainty (e.g. 3.14_)\n' ...
                    '        independent of width of input']))
    end

  case 'displayinfsup'
    setappdata(0,'INTLAB_INTVAL_DISPLAY','DisplayInfsup');
    if see
      disp('===> Default display of intervals by infimum/supremum (e.g. [ 3.14 , 3.15 ])')
    end

  case 'displaymidrad'
    setappdata(0,'INTLAB_INTVAL_DISPLAY','DisplayMidrad');
    if see
      disp('===> Default display of intervals by midpoint/radius (e.g. < 3.14 , 0.01 >)')
    end
    
  %%%%%%%%%%%% get INTLABdata version file %%%%%%%%%%%%
  case 'intlabdata'
    INTLAB_INTLABPATH = getappdata(0,'INTLAB_INTLABPATH');
    version_ = version;
    i = find( version_==' ' );
    if ~isempty(i)
      version_ = version_(1:i-1);
    end
    version_ = [ intvalinit('version') '.' lower(computer) version_ ];
    version_(version_=='.') = '_';
    res = [ INTLAB_INTLABPATH 'INTLAB_Data_' version_ '.mat' ];

  %%%%%%%%%%%%%%%%% get helpfile file %%%%%%%%%%%%%%%%%
  case 'helppfile'
    res = 'helppdata.mat';

  %%%%%%%%%%%%% get real std. fcts. excptn %%%%%%%%%%%%
  case 'realstdfctsexcptn'
    INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
    switch INTLAB_STDFCTS_EXCPTN
      case 0, res = 'RealStdFctsExcptnAuto';
      case 1, res = 'RealStdFctsExcptnWarn';
      case 2, res = 'RealStdFctsExcptnNaN';
      case 3, res = 'RealStdFctsExcptnIgnore';
    end

  %%%%%%%%%%%%% set real std. fcts. excptn %%%%%%%%%%%%
  case 'realstdfctsexcptnauto'
    setappdata(0,'INTLAB_STDFCTS_EXCPTN',0);
    if see
      disp(sprintf([ '===> Complex interval stdfct used automatically for real interval input \n'  ...
                     '         out of range (without warning)']))
    end

  case 'realstdfctsexcptnwarn'
    setappdata(0,'INTLAB_STDFCTS_EXCPTN',1);
    if see
      disp(sprintf([ '===> Complex interval stdfct used automatically for real interval input \n'  ...
             '         out of range, but with warning']))
    end

  case 'realstdfctsexcptnnan'
    setappdata(0,'INTLAB_STDFCTS_EXCPTN',2);
    if see
      disp('===> Result NaN for real interval input out of range ')
    end

  case 'realstdfctsexcptnignore'
    setappdata(0,'INTLAB_STDFCTS_EXCPTN',3);
    if see
      disp('===> !!! Caution: Input arguments out of range are ignored !!!')
      s = sprintf([ '===> Caution: Input arguments out of range are ignored !!! \n' ...
                    '     ===> Use intvalinit(''RealStdFctsExcptnOccurred'') to check whether this happened \n' ...
                    '     ===> Using Brouwer''s Fixed Point Theorem may yield erroneous results (see Readme.txt)']);
      warning(s)
    end

  case 'realstdfctsexcptnoccurred'
    res = getappdata(0,'INTLAB_STDFCTS_EXCPTN_');
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',0);

  %%%%%%%%%%%%%%%% get residual settings %%%%%%%%%%%%%%%%
  case 'residual'
    INTLAB_INTVAL_RESIDUAL = getappdata(0,'INTLAB_INTVAL_RESIDUAL');
    if INTLAB_INTVAL_RESIDUAL==0
      res = 'DoubleResidual';
    elseif INTLAB_INTVAL_RESIDUAL==1
      res = 'ImprovedResidual';
    else
      res = 'QuadrupleResidual';
    end

  %%%%%%%%%%%%%%%% set residual settings %%%%%%%%%%%%%%%%
  case 'doubleresidual'
    setappdata(0,'INTLAB_INTVAL_RESIDUAL',0);
    if see
      disp('===> Double precision residual calculation in verifylss')
    end

  case 'improvedresidual'
    setappdata(0,'INTLAB_INTVAL_RESIDUAL',1);
    if see
      disp('===> Improved residual calculation in verifylss')
    end
    
  case 'quadrupleresidual'
    setappdata(0,'INTLAB_INTVAL_RESIDUAL',2);
    if see
      disp('===> Quadruple precision residual calculation by Dot_ in verifylss')
    end
    
  %%%%%%%%%%%%%%%% get mult. settings %%%%%%%%%%%%%%%%
  case 'ivmult'
    INTLAB_INTVAL_IVMULT = getappdata(0,'INTLAB_INTVAL_IVMULT');
    if INTLAB_INTVAL_IVMULT
      res = 'SharpIVmult';
    else
      res = 'FastIVmult';
    end
    
  %%%%%%%%%%%%%%%% set mult. settings %%%%%%%%%%%%%%%%
  case 'fastivmult'
    setappdata(0,'INTLAB_INTVAL_IVMULT',0);
    if see
      disp(sprintf([ '===> Fast interval matrix multiplication in use (maximum overestimation \n'  ...
                     '        factor 1.5 in radius, see FAQ)']))
    end

  case 'sharpivmult'
    setappdata(0,'INTLAB_INTVAL_IVMULT',1);
    if see
      disp('===> Slow but sharp interval matrix multiplication in use')
    end

  case 'realcomplexassign'
    INTLAB_STDFCTS_RCASSIGN = getappdata(0,'INTLAB_STDFCTS_RCASSIGN');
    switch INTLAB_STDFCTS_RCASSIGN
      case 0, res = 'RealComplexAssignAuto';
      case 1, res = 'RealComplexAssignWarn';
      case 2, res = 'RealComplexAssignError';
    end

  case 'realcomplexassignauto'
    setappdata(0,'INTLAB_STDFCTS_RCASSIGN',0);
    if see
      disp(sprintf([ '===> Real interval arrays automatically transformed to complex array \n'  ...
             '         if a component is assigned a complex value (without warning)']))
    end

  case 'realcomplexassignwarn'
    setappdata(0,'INTLAB_STDFCTS_RCASSIGN',1);
    if see
      disp(sprintf([ '===> Real interval arrays automatically transformed to complex array \n'  ...
             '         if a component is assigned a complex value (with warning)']))
    end

  case 'realcomplexassignerror'
    setappdata(0,'INTLAB_STDFCTS_RCASSIGN',2);
    if see
      disp(sprintf([ '===> An error message is caused when a component of a real interval array \n'  ...
                     '         component is assigned a complex value ']))
    end

  case 'complexinfsupassign'
    INTLAB_INTVAL_CINFSUPASGN = getappdata(0,'INTLAB_INTVAL_CINFSUPASGN');
    if INTLAB_INTVAL_CINFSUPASGN
      res = 'ComplexInfSupAssignWarn';
    else
      res = 'ComplexInfSupAssignNoWarn';
    end
    
  case 'complexinfsupassignwarn'
    setappdata(0,'INTLAB_INTVAL_CINFSUPASGN',1);
    if see
      disp('===> Warning issued when complex interval defined by infsup(zinf,zsup)')
    end

  case 'complexinfsupassignnowarn'
    setappdata(0,'INTLAB_INTVAL_CINFSUPASGN',0);
    if see
      disp('===> No warning issued when complex interval defined by infsup(zinf,zsup)')
    end

  otherwise
    error('intvalinit called with invalid argument')

  end
  
  
  if rndold
    setround(rndold)
  end 
  

function [mthr,IMKL,T] = trysetround(setround_path,param,index)
% Try to add directories dirs and check rounding 
  pathold = path;
  for i=1:length(param)-1
    path( [ setround_path param{i} ] , path );
  end
  % make sure paths are correct
  path(path)
  INTLAB_ROUND_TO_NEAREST = param{end};
  setappdata(0,'INTLAB_ROUND_TO_NEAREST',INTLAB_ROUND_TO_NEAREST);
  [failed,mthr,IMKL] = testrounding(0);
  if failed
    T = inf;
  else
    T = 0; t = 0; % 18.02.2010 VH, new variable t - for appropriate time measures
    imax = 5;
    for i=1:imax
      x = 3.4;
      for j=1:42
        tic
        setround(-1);
        t = toc;
        T = T + toc;
        x = x + 3*5.6;
        tic;
        setround(1);
        t = toc;
        T = T + t;
        x = x - 3*5.6;
        tic;
        setround(0);
        t = toc;
        T = T + t;
      end
    end
    T = T/imax;
  end
      
  path(pathold);		% reset path (no setround)
  
  
function [failed,mthr,IMKL] = testrounding(see) 
%INTVAL extensive test of rounding
%
% failed 0   finished without errors detected
%        1   errors detected: rounding does not work
% IMKL   1   only set if rounding o.k. but IMKL fails
% mthr   1   multi-core rounding works (if multithreading activated)
%        0   errors in multi-core rounding detected
% see    0   suppress error codes
%        1   show errors
%

  failed = 1;
  IMKL = 0;
  mthr = 0;
  
  try     
    % 17.02.2010 VH, first one dummy rounding mode change due to matlab
    % affectation for resetting the first setting of rounding mode
    setround(1);
    getround;
    setround(0);
    
    % various tests of the rounding function
    setround(0)
    clear functions
    e=eps*eps;
    if (1+e)~=1
      if see, disp('error in rounding to nearest 1'), end
      return
    end
    if (-1+e)~=-1
      if see, disp('error in rounding to nearest 2'), end
    end

    setround(1)
    clear functions
    if (1+e)<=1
      if see, disp('error in rounding to up 1    '), end
      return
    end
    if (-1+e)<=-1
      if see, disp('error in rounding to up 2    '), end
      return
    end

    setround(-1)
    clear functions
    if (1-e)>=1
      if see, disp('error in rounding to down 1 '), end
      return
    end
    if (-1-e)>=-1
      if see, disp('error in rounding to down 2 '), end
      return
    end

    setround(2)
    clear functions
    if (1-e)>=1
      if see, disp('error in rounding to zero 1 '), end
      return
    end
    if (-1+e)<=-1
      if see, disp('error in rounding to zero 2 '), end
      return
    end

%    setround(2)
%    clear functions
%    eval('6*7')
%    if (1-e)>=1
%      if see, disp('error in rounding to zero 1.1 '), end
%      return
%    end
%    if (-1+e)<=-1
%      if see, disp('error in rounding to zero 2.1 '), end
%      return
%    end

    setround(0)
    clear functions
    if (1+e)~=1
      if see, disp('error in rounding to nearest 3'), end
      return
    end
    if (-1+e)~=-1
      if see, disp('error in rounding to nearest 4'), end
      return
    end

    % extra test for BLAS library
    x = 0.1*ones(2,1);
    setround(-1)
    clear functions
    dotinf = x'*x;
    setround(1)
    clear functions
    dotsup = x'*x;
    setround(0)
    clear functions
    if dotinf==dotsup
      IMKL = 1;
      if see, disp('error in dot product 5'), end
      return
    end
    
    % tests succeeded
    mthr = 1;
    
    % extra test for multi-core/multi-threaded rounding
    try 
      wngThread = warning('off','MATLAB:maxNumCompThreads:Deprecated');
      nthr = maxNumCompThreads;   % not available in later Matlab versions
      warning(wngThread.state,'MATLAB:maxNumCompThreads:Deprecated');
    catch 
        nthr = 2; 
        lasterr('');
    end
    if nthr > 1        
        n = 220; 
        a = rand(n);
        setround(1);
        clear functions
        c1 = a*a;
        setround(-1);
        clear functions
        c2 = a*a;
        if any(c1(:) == c2(:))
            mthr = 0;
        end
    end
    
    % check change of rounding mode
    for i=1:5
      for n=[1 5 30 100]
        a = randn(n);
        g = randint(4)-2;
        setround(g);
        clear functions
        a*intval(a);
        midrad(a,0.1)*intval(a);
        gg = getround;
        if g~=gg
          mthr = 0;
          if see, disp('error in changing rounding back and forth'), end
          return
        end
      end
    end
    
    failed = 0;     % Test succeeded
            
  catch
    mthr = 0;
    lasterr('');
  end
  
  
function INTLAB_INTVAL_POWER10 = initpower10
%Initialization of INTLAB_INTVAL_POWER10
%
%For integer e in [-340,309] and integer m in [1,9] it is
%
%  INTLAB_INTVAL_POWER10.inf(m,E)  <= m*10^e <=  INTLAB_INTVAL_POWER10.sup(m,E)
%
%where E:=e+341 and binary numbers INTLAB_INTVAL_POWER10.inf and
%  INTLAB_INTVAL_POWER10.sup are best possible.
%For fast computations, INTLAB_INTVAL_POWER10.sup(1) := 0 !
%
%Implementation uses a simplified long arithmetic with long representing
%
%    sum( i , long(i) * base^i )
%
%for base = 2^47 and 0 <= long(i) < 2*base
%

  setround(0)

  baseExponent = 47;
  base = 2^baseExponent;
  imax = 4;
  emin = -340;
  emax = 309;
  INTLAB_INTVAL_POWER10.inf = zeros(9,emax-emin+1);
  INTLAB_INTVAL_POWER10.sup = zeros(9,emax-emin+1);

  long = zeros(1,23+imax);   % exponents >= 0
  long(23) = 1;
  first = 23;
  for e=0:emax
    for m=1:9

      long_ = m*long;

      setround(-1)
      s = 0;
      for i=imax:-1:0
        s = s/base + long_(first+i) ;
      end
      INTLAB_INTVAL_POWER10.inf(m,e-emin+1) = ...
        s * 2^( baseExponent*(23-first) );
      setround(1)
      s = any(long_(first+imax+1:end)~=0)/base;
      for i=imax:-1:0
        s = s/base + long_(first+i) ;
      end
      INTLAB_INTVAL_POWER10.sup(m,e-emin+1) = ...
        s * 2^( baseExponent*(23-first) );
      setround(0)
    end

    long = 10*long;
    q = floor(long/base);
    long = long - q*base;
    long(1:end-1) = long(1:end-1) + q(2:end);
    first = first - (long(first-1)~=0);
  end

  lmax = 30;

  long = zeros(1,lmax);   % exponents < 0, rounding downwards
  long(1) = 1;
  first = 1;
  for e=-1:-1:emin

    for i=first:lmax-1
      q = floor(long(i)/10);
      long(i+1) = long(i+1) + (long(i)-10*q)*base;
      long(i) = q;
    end
    long(lmax) = floor(long(lmax)/10);
    first = first + (long(first)==0);

    for m=1:9

      long_ = m*long;

      setround(-1)
      s = 0;
      for i=imax:-1:0
        s = s/base + long_(first+i) ;
      end
      INTLAB_INTVAL_POWER10.inf(m,e-emin+1) = ...        % avoid underflow
        s * 2^( baseExponent*(10-first) ) * 2^( baseExponent*(-9) );
      setround(0)
    end
  end

  long = zeros(1,lmax);   % exponents < 0, rounding upwards
  long(1) = 1;
  first = 1;
  for e=-1:-1:emin

    for i=first:lmax-1
      q = floor(long(i)/10);
      long(i+1) = long(i+1) + (long(i)-10*q)*base;
      long(i) = q;
    end
    long(lmax) = floor(long(lmax)/10) + 1;
    first = first + (long(first)==0);

    for m=1:9

      long_ = m*long;

      setround(1)
      s = 1/base;
      for i=imax:-1:0
        s = s/base + long_(first+i) ;
      end
      INTLAB_INTVAL_POWER10.sup(m,e-emin+1) = ...        % avoid underflow
        s * 2^( baseExponent*(10-first) ) * 2^( baseExponent*(-9) );
      setround(0)
    end
  end

  INTLAB_INTVAL_POWER10.sup(1) = 0;
  INTLAB_INTVAL_POWER10.inf(3056) = 0.5;
  INTLAB_INTVAL_POWER10.sup(3056) = 0.5;
  
  
function versioncheck(ver)
% check for known errors in various Matlab versions
  Ver = version;
  Ver = Ver(end-6:end-1);
  if isempty(ver)
    % known bug in R2006a
    n = 10000; A = sprand(n,n,2/n/n);
    tic, 2*A; t=toc;
    if t>0.3
      intvalmessages('R2006a',Ver)
      % prevent Matlab from hard stop when using SSE2 rounding routines
      setenv('KMP_DUPLICATE_LIB_OK','TRUE');
    end
    % known bug in R2006b
    x=sparse([1 -1]); y=sparse([0 0]); z=min(x,y);
    if all(z>-0.5)
      intvalmessages('R2006b',Ver)
      % prevent Matlab from hard stop when using SSE2 rounding routines
      setenv('KMP_DUPLICATE_LIB_OK','TRUE');
    end
    % known bug in R2009b
    A=[ 0 2 ; 1 0 ]; b=[ 2 ; 0 ]; x=A'\b;
    if norm(x-[0;2])>1e-4
      intvalmessages('R2009b',Ver)
      % prevent Matlab from hard stop when using SSE2 rounding routines
      setenv('KMP_DUPLICATE_LIB_OK','TRUE');
    end
  end
  if isequal(Ver,ver)
    switch ver
      case 'R2008b'
        intvalmessages(ver)
        % prevent Matlab from hard stop when using SSE2 rounding routines
        setenv('KMP_DUPLICATE_LIB_OK','TRUE');
      case 'R2009a'
        intvalmessages(ver)
        % prevent Matlab from hard stop when using SSE2 rounding routines
        setenv('KMP_DUPLICATE_LIB_OK','TRUE');
    end
  end

  
function intvalmessages(warncase,param)
%INTVALMESSAGES   some warn and error messages of INTLAB
%
%   intvalinit(param,see)
%
%possible values for warncase:
%
%  'FindBestRounding'        Warning for possible hard stop
%
%  'BestRoundingTestFailed'  Hard stop occurred
%
%  'IMKLerror'               Prints error message if errors in switching 
%                            rounding mode in BLAS detected
%
%  'ROUNDerror'              Prints error message if errors in switching
%                            rounding mode detected
%
%  'FatalError'              No way to establish stable switching of rounding mode
%
%  'R2006a'                  Prints warning message when using Matlab 2006a
%
%  'R2006b'                  Prints error message when using Matlab 2006b
%
%  'R2008b'                  Prints error message when using Matlab 2008b
%
%  'R2009a'                  Prints error message when using Matlab 2009a
%
%  'R2009b'                  Prints error message when using Matlab 2009b
%
%  'MULTICOREerror'          Problems with switching rounding mode on multi-core computers
%
%  'LIC'                     License information
%
%  'Reference'               Reference information (first publication)
%
%  'WELCOME'                 First "Welcome" message
%
%  'CONSTerror'              Initialization of constants for verified standard functions failed
%
%  'STDFCNmissed'            Invalid call of INTLAB (constants for verified standard functions
%                            missing)
%
switch warncase
  case 'FindBestRounding'
    disp(' ')
    disp('VVVVVVVVVVVVVVVV Installation of INTLAB VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV')
    disp(' ')
    disp('This is the first call of INTLAB, and the best way to switch the ')
    disp('rounding mode will be detected. It is a very rigorous check and takes')
    disp('a short moment. This is done only once with the very first call of INTLAB.')
    disp(' ')
    disp('To find the best routine, various routines are checked. ')
    disp('Unfortunately, in very rare cases it might happen that a hard error')
    disp('occurs. This should never happen, but it happes in some Matlab versions')
    disp('and/or operating systems, and it is beyond my control. ')
    disp(' ')
    disp('In such a rare case of a hard error, please restart INTLAB.')
    disp('The circumstances will be stored, and this switching routine will not')
    disp('be called again. Sorry for possible inconveniences. ')
    disp(' ')
    input('Press Enter to continue.')
    disp(' ')
  case 'BestRoundingTestFailed'
    disp(' ')
    disp('$$$$$$$$$$$$$$$$ Installation of INTLAB $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    disp(' ')
    disp('Unfortunately the rare case of a hard stop happened using')
    disp([ blanks(5) param{:} ])
    disp('This switch of rounding mode is excluded from further tries. ')
    disp(' ')
    disp('The fastest routine for switching the rounding mode is to be detected.')
    disp('A hard error might be encountered again; in that case please restart ')
    disp('INTLAB afresh. The best routine will be stored on the disc, and ')
    disp('this embarrassing procedure will not happen again. ')
    disp(' ')
    disp('Sorry for the inconvenience, this is beyond my control. ')
    disp(' ')
    input('Press Enter to continue.')
    disp(' ')
  case 'IMKLerror'
    disp('*********************************************************************')
    disp('*************** Installation of INTLAB ******************************')
    disp('**                                                                 ')
    disp('** Errors in switching of rounding mode in BLAS detected!          ')
    disp('** Any use of interval routines may produce **erreneous** results! ')
    disp('** This happens for some environments under Matlab R14, service    ')
    disp('** 1, 2 and 3, R2008b and R2009a.                                  ')
    disp('**                                                                 ')
    disp('** Please try to change the environment to use the ATLAS library by')
    disp('**   setting the corresponding environment variable BLAS_VERSION.  ')
    disp('** Under Windows: In Start -> Control Panel -> System ->           ')
    disp('**                       Advanced -> Environment Variables         ')
    disp('**    click the new button for "System Variables" and define a     ')
    disp('**    new environment variable BLAS_VERSION to be                  ')
    disp('**       atlas_Athlon.dll, atlas_P4.dll, atlas_PIII.dll ...        ')
    disp('**    where the appropriate atlas library for your processor is in ')
    disp('**    <Matlab>\bin\win32 or can be downloaded from the internet.   ')
    disp('**    As reported by Bastian Ebeling you may as well define the    ')
    disp('**    new environment variable BLAS_VERSION to be                  ')
    disp('**       refblas.dll                                               ')
    disp('**    which works in newer releases of Matlab where the ATLAS      ')
    disp('**    library is no longer included. It is, however, pretty slow.  ')
    disp('** See file FAQ.txt for more information.                          ')
    disp('**                                                                 ')
    disp('** Sorry for the inconvenience, but there is no way to do it       ')
    disp('**   within Matlab.                                                ')
    disp('**                                                                 ')
    disp('** Press Enter to terminate INTLAB                                 ')
    disp('**                                                                 ')
  case 'ROUNDerror'
    disp('*********************************************************************')
    disp('************** Installation of INTLAB *******************************')
    disp('**                                                                 ')
    disp('** Errors in switching of rounding mode detected!                  ')
    disp('**                                                                 ')
    disp('** Any use of interval routines may produce **erreneous** results! ')
    disp('**                                                                 ')
    disp('** For Matlab 5.3, R11 or older under Windows:                     ')
    disp('**   Probably the .dll file for function "setround" is not         ')
    disp('**   installed or in the wrong path.                               ')
    disp('**                                                                 ')
    disp('** For other operating systems, provide file for switching         ')
    disp('**   rounding mode (see INTLAB homepage).                          ')
    disp('**                                                                 ')
    disp('** Press Enter to terminate INTLAB                                 ')
    disp('**                                                                 ')
    disp(' ')
    disp(' ')
  case 'FatalError'
    disp('*********************************************************************')
    disp('************ Installation of INTLAB *********************************')
    disp('**                                                                 ')
    disp('** Errors in switching of rounding mode detected! Changing the     ')
    disp('** rounding mode is not stable. This is beyond my control, please  ')
    disp('** use another Matlab version.                                     ')
    disp('**                                                                 ')
    disp('** Any use of interval routines may produce **erreneous** results! ')
    disp('**                                                                 ')
    disp('** Press Enter to terminate INTLAB                                 ')
    disp('**                                                                 ')
    disp(' ')
    disp(' ')
  case 'R2006a'
    disp('************************************************************************')
    disp('*!*!*!*!*!*!* Installation of INTLAB *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    disp('*!                                                                       ')
    disp(['*! You are using Matlab ' param ', which has a flaw in sparse matrices: '])
    disp('*! The commands                                                          ')
    disp('*!   n = 20000; A = sprand(n,n,2/n/n), tic, 2*A, toc                     ')
    disp('*! take unexpectedly long and deliver something like                     ')
    disp('*!   A =                                                                 ')
    disp('*!         (7058,198)           0.2028                                   ')
    disp('*!        (16264,2778)          0.1987                                   ')
    disp('*!   ans =                                                               ')
    disp('*!         (7058,198)           0.4055                                   ')
    disp('*!        (16264,2778)          0.3974                                   ')
    disp('*!   Elapsed time is 3.824998 seconds.                                   ')
    disp('*!                                                                       ')
    disp(['*! So consider whether you want to use ' param '.                      '])
    disp('*!                                                                       ')
    disp('*! Press Enter to continue.                                              ')
    disp('*!                                                                       ')
    disp(' ')
    input(' ');
    disp(' ')
  case 'R2006b'
    disp('************************************************************************')
    disp('*!*!*!*!*!*!* Installation of INTLAB *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    disp('*!                                                                       ')
    disp(['*! You are using Matlab ' param ', which has a severe bug in sparse matrices:'])
    disp('*! The commands                                                          ')
    disp('*!   x = sparse([1 -1]), y = sparse([0 0]), min(x,y)                     ')
    disp('*! deliver the wrong result                                              ')
    disp('*!   x =                                                                 ')
    disp('*!      (1,1)        1                                                   ')
    disp('*!      (1,2)       -1                                                   ')
    disp('*!   y =                                                                 ')
    disp('*!      All zero sparse: 1-by-2                                          ')
    disp('*!   ans =                                                               ')
    disp('*!      All zero sparse: 1-by-2                                          ')
    disp('*!                                                                       ')
    res = input('If you accept incorrect results, input y, otherwise Matlab will terminate ','s');
    if isequal(lower(res),'y')
      disp('*!                                                                     ')
      disp('*! WARNING: results by INTLAB may be incorrect.                        ')
      disp('*!                                                                     ')
      disp(' ')
      disp(' ')
    else
      exit
    end
  case 'R2008b'
    disp('************************************************************************')
    disp('*!*!*!*!*!*!*!* Installation of INTLAB *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    disp('*!                                                                       ')
    disp('*! You are using Matlab R2008b, which has - on some operating systems -  ')
    disp('*! severe problems in changing the rounding mode.                        ')
    disp('*! For example, a simple output seems to reset the control word to some  ')
    disp('*! default, jeopardizing the rounding switch.                            ')
    disp('*! Changing the rounding mode will be checked next, but it may not work. ')
    disp('*! In that case, if possible, please use another Matlab release.         ')
    disp('*!                                                                       ')
    disp('*! Press Enter to continue.                                              ')
    disp('*!                                                                       ')
    disp(' ')
    input(' ');
    disp(' ')
  case 'R2009a'
    disp('************************************************************************')
    disp('*!*!*!*!*!*!*!*!* Installation of INTLAB *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    disp('*!                                                                       ')
    disp('*! You are using Matlab R2009a, which has - on some operating systems -  ')
    disp('*! severe problems in changing the rounding mode.                        ')
    disp('*! For example, a simple output seems to reset the control word to some  ')
    disp('*! default, jeopardizing the rounding switch.                            ')
    disp('*! Changing the rounding mode will be checked next, but it may not work. ')
    disp('*! In that case, if possible, please use another Matlab release.         ')
    disp('*!                                                                       ')
    disp('*! Press Enter to continue.                                              ')
    disp('*!                                                                       ')
    disp(' ')
    input(' ');
    disp(' ')
  case 'R2009b'
    disp('************************************************************************')
    disp('*!*!*!*!*!*!*!* Installation of INTLAB *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    disp('*!                                                                       ')
    disp(['*! You are using Matlab ' param ', which in turn uses Intel® MKL for matrix '])
    disp('*! and vector operation. For multi-core computers there are severe bugs  ')
    disp('*! in ordinary rounding to nearest, see                                  ')
    disp('*!   http://software.intel.com/en-us/articles/dgemm-and-sgemm-accuracy/  ')
    disp('*! Frank Schmidt from TU Chemnitz reported the following problem:        ')
    disp('*!    A=[ 0 2 ; 1 0 ]; b=[ 2 ; 0 ]; x=A''\b                              ')
    disp('*!    x =                                                                ')
    disp('*!         0                                                             ')
    disp('*!         1                                                             ')
    disp('*! The true solution is obviously [0;2].                                 ')
    disp('*!                                                                       ')
    res = input('If you accept incorrect results, input y, otherwise Matlab will terminate ','s');
    if isequal(lower(res),'y')
      disp('*!                                                                     ')
      disp('*! WARNING: results by INTLAB may be incorrect.                        ')
      disp('*!                                                                     ')
      disp(' ')
      disp(' ')
    else
      exit
    end
  case 'MULTICOREerror'
    disp('*********************************************************************')
    disp('************* Installation of INTLAB ********************************')
    disp('**                                                                 ')
    disp('** Errors in switching of rounding mode detected!                  ')
    disp('**                                                                 ')
    disp('** Any use of interval routines may produce **erreneous** results! ')
    disp('**                                                                 ')
    disp('** If you are using Matlab R2008a or earlier and a multicore       ')
    disp('**   machine, please set the number of computational threads to 1  ')
    disp('**   in the "file -> preference menue".                            ')
    disp('**   Otherwise, add "-singleCompThread" as a startup option.       ')
    disp('**                                                                 ')
    disp('** Press Enter to terminate INTLAB                                 ')
    disp('**                                                                 ')
    disp(' ')
    disp(' ')
  case 'LIC'
    disp(' ')
    disp('=========================================================================')
    disp('*****  Please provide proper reference to                                ')
    disp('*****     S.M. Rump: INTLAB - INTerval LABoratory. In Tibor Csendes,     ')
    disp('*****       editor, Developments in Reliable Computing, pages 77-104.    ')
    disp('*****       Kluwer Academic Publishers, Dordrecht, 1999.                 ')
    disp('*****                                                                    ')
    disp('*****  Commercial use or use in conjunction with a commercial program    ')
    disp('*****  which requires INTLAB or part of INTLAB to function properly      ')
    disp('*****  requires a special license. Please contact me  rump (at) tuhh.de  ')
    disp(' ')
    disp(' ')
  case 'reference'
    disp(' ')
    disp('+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+')
    disp('+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+')
    disp(':: ')
    disp('::  Some references to papers using INTLAB are collected on the     ')
    disp('::    INTLAB homepage   http://www.ti3.tuhh.de/                     ')
    disp('::  If you have additional references to add, please send me        ')
    disp('::  a mail ( rump [at] tuhh.de )                                    ')
    disp('::                                                                  ')
    disp('::  In any publication or other material using INTLAB please cite   ')
    disp('::     S.M. Rump: INTLAB - INTerval LABoratory. In Tibor Csendes,   ')
    disp('::       editor, Developments in Reliable Computing, pages 77-104.  ')
    disp('::       Kluwer Academic Publishers, Dordrecht, 1999.               ')
    disp('::  The bibtex entry is                                             ')
    disp('    @incollection{Ru99a,                                            ')
    disp('       author = {Rump, {S.M.}},                                     ')
    disp('       title = {{INTLAB - INTerval LABoratory}},                    ')
    disp('       editor = {Tibor Csendes},                                    ')
    disp('       booktitle = {{Developments~in~Reliable Computing}},          ')
    disp('       publisher = {Kluwer Academic Publishers},                    ')
    disp('       address = {Dordrecht},                                       ')
    disp('       pages = {77--104},                                           ')
    disp('       year = {1999},                                               ')
    disp('       note = {\url{http://www.ti3.tuhh.de/rump/}}                  ')
    disp('    }                                                               ')
    disp('::                                                                  ')
    disp(' ')
    disp(' ')
  case 'WELCOME'
    disp(' ')
    disp('**************************************************************')
    disp(['*** Welcome to INTLAB - INTerval LABoratory Version ' param ])
    disp('***   The Matlab toolbox for Reliable Computing            ')
    disp('***   Siegfried M. Rump, Institute for Reliable Computing  ')
    disp('***   Hamburg University of Technology, Germany            ')
    disp(' ')
    disp(' ')
    disp('INTLAB uses several control variables; to check the default values call "intlabsetting"')
  case 'CONSTerror'
    disp('*********************************************************************')
    disp('************** Installation of INTLAB *******************************')
    disp('**                                                                 ')
    disp('** *** Error:                                                      ')
    disp('** Constants for rigorous standard functions could not be computed ')
    disp('**                                                                 ')
    disp('** This is extremely rare; within years and several thousand users ')
    disp('**   it happened only once.                                        ')
    disp('**                                                                 ')
    disp('** Please report this to                                           ')
    disp('**                                                                 ')
    disp('**   rump (at) tuhh.de                                             ')
    disp('**                                                                 ')
    disp('** Thanks for cooperation.                                         ')
    disp('**                                                                 ')
    disp('** Press Enter to terminate INTLAB.                                ')
    disp('**                                                                 ')
    disp(' ')
    disp(' ')
  case 'STDFCNmissed'
    disp(' ')
    disp('**** Data for rigorous INTLAB standard functions have to be generated, ')
    disp('**** takes a minute or so.                                             ')
    disp('**** The data will be generated _once_ during installation and stored. ')
    disp(' ')
    input('Press Enter for unique generation of data. ')
    disp(' ')
end

