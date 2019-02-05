function startintlab
%STARTINTLAB  Initialization of paths, INTLAB constants, etc.
%             Adapt this to your local needs
%

% written  10/16/98   S.M. Rump
% modified 11/05/98   allow blanks in path name, Max Jerrell, Flagstaff
% modified 10/23/99   S.M. Rump  change directory before intvalinit('init'),
%                                work directory
% modified 12/15/01   S.M. Rump  displaywidth added
% modified 12/08/02   S.M. Rump  reset directory to old dir (Matlab problems under unix)
% modified 12/18/02   S.M. Rump  Hessians added
% modified 04/04/04   S.M. Rump  working directory added, clear removed, changed to function
%                                  adapted to new functionalities
% modified 11/09/05   S.M. Rump  working dir and cd(INTLABPATH) removed, INTLABPATH computed
% modified 02/11/06   S.M. Rump  SparseInfNanFlag removed
% modified 06/15/06   S.M. Rump  directory name corrected (thanks to Jaap van de Griend)
% modified 05/09/07   S.M. Rump  path corrected (thanks to Bastian Ebeling)
% modified 09/10/07   S.M. Rump  INTLAB path
% modified 11/26/07   S.M. Rump  ordering of default settings
% modified 02/19/08   S.M. Rump  rename setround paths, avoid change of code
% modified 05/19/08   S.M. Rump  reorganization of setround
% modified 05/26/08   S.M. Rump  rehash added
% modified 05/09/09   S.M. Rump  directory AccSumDot added, warning for "format" supressed
% modified 05/25/09   S.M. Rump  directory taylor added
% modified 01/17/10   S.M. Rump  autoamd switched off
% modified 04/08/10   S.M. Rump  extra checkrounding omitted
% modified 08/22/12   S.M. Rump  default of stdfcts NaN, warnings may be switched off
% modified 08/24/12   S.M. Rump  STAGE superflous due to redesign of verifylss
% modified 10/01/12   S.M. Rump  comments
% modified 10/03/12   S.M. Rump  output at start
% modified 10/03/12   S.M. Rump  SparseArrayDeriv removed for gradient, hessian and slope
% modified 10/13/12   S.M. Rump  demos.html added
% modified 11/06/12   S.M. Rump  warning for format switched off
% modified 12/07/12   S.M. Rump  filesep to add html path
%

  try											      % take care of older Matlab versions
    warning('off','MATLAB:dispatcher:nameConflict')   % take care of 'format'
  end
%   wng = warning;  
%   warning off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default setting:  please change these lines to your needs %%
%VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

% Intlab directory

% Defines INTLABPATH to be the directory, in which this file "startintlab" is contained
  dir_startintlab = which('startintlab');
  INTLABPATH = dir_startintlab(1:end-13);

% If INTLAB is contained in another directory, please uncomment and adapt this line
% INTLABPATH = 'C:\INTLAB VERSIONS\INTLAB\';           % blanks allowed 


%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% add new paths
  addpath( [ INTLABPATH ]                       , ...
           [ INTLABPATH 'intval' ]              , ...
           [ INTLABPATH 'gradient' ]            , ...
           [ INTLABPATH 'hessian' ]             , ...
           [ INTLABPATH 'taylor' ]              , ...
           [ INTLABPATH 'slope' ]               , ...
           [ INTLABPATH 'polynom' ]             , ...
           [ INTLABPATH 'utility' ]             , ...
           [ INTLABPATH 'long' ]                , ...
           [ INTLABPATH 'accsumdot' ]           , ...
           [ INTLABPATH 'demos' ]               , ...
           [ INTLABPATH 'demos' filesep 'html' ])
  path(path)			         % make sure paths are correct

%%%%%%%%%% set INTLAB environment

  format compact
  format short
  try
    feature jit off                % avoid change of code
  end
  try
    feature accel off		    
  end


%%%%%%%%%% initialize sparse systems

  spparms('autoamd',0);            % switch off automatic band reduction
  spparms('autommd',0);            % ( switching on may slow down sparse
                                   %   computations significantly )


%%%%%%%%%% initialize interval toolbox (see "help intvalinit")

  intvalinit('Init',0,INTLABPATH)  % Initialize INTLAB constants, see=0


%%%%%%%%%% initialize gradient toolbox

  gradientinit


%%%%%%%%%% initialize Hessian toolbox

  hessianinit


%%%%%%%%%% initialize slope toolbox

  slopeinit


%%%%%%%%%% initialize long toolbox

  longinit('Init')


%%%%%%%%%% initialize polynom toolbox

polynominit('Init')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Default setting:  change these lines to your needs   %%%%%%%%%%%%%%
%VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

  see = 0;
  displaywidth(120,see);               % width of output

  intvalinit('Display_',see)
  intvalinit('RealStdFctsExcptnNaN',see)
  intvalinit('ImprovedResidual',see)
  intvalinit('FastIVMult',see)
  intvalinit('RealComplexAssignAuto',see)
  intvalinit('ComplexInfSupAssignWarn',see)

  longinit('WithErrorTerm',see)

  polynominit('DisplayUPolySparse',see)
  polynominit('EvaluateUPolyHorner',see)
  polynominit('EvaluateMPolyPower',see)
  polynominit('AccessVariableWarn',see)

  sparsegradient(50,see);     % Gradients stored sparse for fifty and more variables

  sparsehessian(10,see);      % Hessians stored sparse for ten and more variables

%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  intvalinit('license')
  intlablogo([0:4:32 32:-4:24 24:4:32])
  pause(0.5)
  close
  
  
%%%%%%%%%% store current setting as default setting
  
  intvalinit('StoreDefaultSetting')
  
  
%%%%%%%%%% set working environment

  setround(0)             % set rounding to nearest
%   warning(wng)

% uncomment and adapt this statement if necessary
% cd('c:\rump\matlab\work')
  