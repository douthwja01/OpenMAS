function setround(rnd)
%SETROUND     Switch rounding mode
%
%   setround(rnd)
%
%For  rnd = -1  switch rounding downwards (towards -inf)
%     rnd =  0  switch rounding to nearest
%     rnd =  1  switch rounding upwards (towards inf)
%     rnd =  2  switch rounding towards zero (chop)
%
%INTLAB switches automatically between different ways to change the
%  rounding mode. On some operating systems and versions of Matlab
%  everything is done in Matlab using the Matlab-routine 'feature'
%  or 'system_dependent', or some assembly file is used. 
%For some architectures assembly files can be found on our home page.
%The switch should work automatically.
%
%It is assumed that the processor is permanently switched into the
%specified rounding mode, i.e. all subsequent operations are performed
%according to this rounding mode. When invoking INTLAB, always
%
%  intvalinit('CheckRounding')
%
%is called to be sure the preceding statement is true. Apparently, there
%are difficulties with DEC Alpha workstations in which op-codes carry
%an individual rounding mode.
%

% written  11/30/98     S.M. Rump
% modified 12/18/99     S.M. Rump  system dependent call for Matlab V5.3.1f
%                                  under Windows added
% modified 02/15/02     S.M. Rump  rounding to nearest changed to 0.5 (should work under Linux, too, 
%                                     thanks to Dr. Jaap A. van de Griend)
% modified 10/09/02     S.M. Rump  Rounding switch by global variable
% modified 11/02/05     S.M. Rump  Improved performance (thanks to Jörg Kubitz, Hannover)
% modified 05/05/08	    T. Ogita   Tentative solution for using Intel MKL BLAS
% modified 05/19/08	    S.M. Rump  Reorganization of setround
% modified 08/26/12     S.M. Rump  global variables removed
%

% The following should work under Matlab 6+, also under Linux operating system.
%

  INTLAB_ROUND_TO_NEAREST = getappdata(0,'INTLAB_ROUND_TO_NEAREST');
  
  switch rnd
    case  0, system_dependent('setround',INTLAB_ROUND_TO_NEAREST)         % round to nearest
    case -1, system_dependent('setround',-inf)        % round towards -inf
    case  1, system_dependent('setround',inf)         % round towards +inf
    case  2, system_dependent('setround',0)           % round towards zero (chop)
  otherwise
    error('invalid input argument to setround')
  end
