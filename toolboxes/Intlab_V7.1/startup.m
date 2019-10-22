%STARTUP      Dummy routine to call startintlab 
%
%Insert the call to "startintlab" into your personal startup-file.
%

% written  09/10/00     S.M. Rump  
% modified 05/19/08     S.M. Rump  catch error from BLAS_VERSION
% modified 06/18/10     S.M. Rump  catch wrong directory, thanks to Kaori Nagatou
% modified 09/07/12     S.M. Rump  comment
% modified 06/24/13     S.M. Rump  comment
%
currentDir = pwd;

  cont = 1;
  try
    % Define the IntLab directory path
    intlabPath = mfilename('fullpath');
    fileName = strcat('\',mfilename);
    intDir = erase(intlabPath,fileName);
    % Move to the path for intialisation
    fprintf('Moving to interval directory: %s\n',intDir);    
    cd(intDir)
  catch
    cont = 0;
    disp('=========================================================================')
    disp('=========================================================================')
    disp('*** Please specify in "startup.m" in the line                         ***')
    disp('***     cd ''c:\intlab''                                                ***')
    disp('*** the directory where INTLAB is installed.                          ***')
    disp('*** Otherwise you might have old  .mat-files in the INTLAB directory. ***')
    disp('*** If so, please remove them.     								    ***')
    disp('***                                                                   ***')
    disp('*** Press Enter to continue, and edit "startup.m".                    ***')
    disp('*** !!! INTLAB is NOT yet working properly !!!                        ***')
    disp('=========================================================================')
    disp('=========================================================================')
    pause
  end
  
  if cont
    try
      startintlab
    catch
      disp('========================================================================')
      disp('========================================================================')
      disp('*** Some error occurred during startup, this should not happen.      ***')
      disp('*** A possible reason is that startintlab.m is not in your           ***')
      disp('***   current directory.                                             ***')
      disp('***                                                                  ***')
      disp('*** Another reason might be that you used INTLAB before and set the  ***')
      disp('***   environment variable BLAS_VERSION.                             ***')
      disp('*** In previous versions of INTLAB this was necessary because        ***')
      disp('***   Intel Math Kernel Library (IMKL) did not work with my          ***')
      disp('***   switching of the rounding mode.                                ***')
      disp('*** Meanwhile this is changed, IMKL should work properly. So please  ***')
      disp('***   delete the environment variable BLAS_VERSION and try again.    ***')
      disp('*** Sorry for inconvenience.                                         ***')
      disp('***                                                                  ***')
      disp('*** If nothing helps, please go step by step through startup and     ***')
      disp('***   startintlab to identify the reason.                            ***')
      disp('========================================================================')
      disp('========================================================================')
    end
  end
cd(currentDir)