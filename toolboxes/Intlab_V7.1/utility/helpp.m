function helpp(name,all,level)
%HELPP        Implementation of intelligent help
%
%For file being a string,
%    helpp(file)               searches file amongst INTLAB functions
%    helpp(file,'user')        searches file in user-defined path
%    helpp(file,'all')         searches all accessible files
%    helpp(file,level)         search with specified level (default 6) amongst INTLAB functions
%    helpp(file,'user',level)  search with specified level (default 6) in user-defined path
%    helpp(file,'all',level)   search with specified level (default 6) in all accessible files
%    helpp                     build matrix memory according to Matlab path
%
%First,
%   helpp 
%without parameters must be called for initialization. 
%
%Then, for example,
%   helpp intersection
%produces in my environment
%  ans =
%  level  7  c:\intlab versions\intlab\intlabsetting.m                      
%  level 10  c:\intlab versions\intlab\intval\@intval\intersect.m     
%
%or 
%
%  helpp tschebischeff
%produces
%  ans =
%  level 6  c:\intlab versions\intlab\utility\fletcher.m                   
%  level 9  c:\intlab versions\intlab\polynom\chebyshev.m                  
%  level 9  c:\intlab versions\intlab\polynom\chebyshev2.m  
%
%Intelligent help does work good on short names, the relative distance becomes to small.
%If too much (or not enough) output is produced, the level may be increased (or decreased).
%

% written  07/08/03     S.M. Rump
% modified 06/10/05     S.M. Rump  user-defined path corrected
% modified 05/12/06     S.M. Rump  file changed
% modified 07/01/07     S.M. Rump  USER file saved, minlevel
% modified 08/26/12     S.M. Rump  global variables removed
%

  if nargin==0
    
    setappdata(0,'INTLAB_MATRIX_MEMORY_ALL',sparse([],[],[]));              % matrix memory for complete path (to be built)
    setappdata(0,'INTLAB_MATRIX_FILES_ALL',[]);
    setappdata(0,'INTLAB_MATRIX_MEMORY_USER',sparse([],[],[]));             % matrix memory for user-defined path (to be built)
    setappdata(0,'INTLAB_MATRIX_FILES_USER',[]);
    setappdata(0,'INTLAB_MATRIX_MEMORY',sparse([],[],[]));                  % matrix memory for user-defined path (to be built)
    setappdata(0,'INTLAB_MATRIX_FILES',[]);
    
    dirold = cd;
    P = [path pathsep];
    
    p = which('spones');            % compute Matlab path
    i = findstr(p,['toolbox' filesep]);
    MatlabPath = p(1:i-1);
    
    p = which('intvalinit');        % compute INTLAB path
    i = findstr(p,['intval' filesep]);
    INTLABPath = p(1:i-1);
    
    while ~isempty(P)
      
      i = findstr(P,pathsep);
      p = P(1:i-1);                 % search directory
      
      % add information of path p to matrix memory for complete path
      disp(['adding path ' p ' to matrix memory for complete path'])
      [INTLAB_MATRIX_MEMORY_ALL,INTLAB_MATRIX_FILES_ALL] = ...
        append_matrix_memory(p,filesep,INTLAB_MATRIX_MEMORY_ALL,INTLAB_MATRIX_FILES_ALL);
      setappdata(0,'INTLAB_MATRIX_MEMORY_ALL',INTLAB_MATRIX_MEMORY_ALL);
      setappdata(0,'INTLAB_MATRIX_FILES_ALL',INTLAB_MATRIX_FILES_ALL);
      
      % check whether path p belongs not to Matlab path
      if isempty(findstr(p,MatlabPath))
        
        disp(['adding path ' p ' to matrix memory for user-defined path'])
        [INTLAB_MATRIX_MEMORY_USER,INTLAB_MATRIX_FILES_USER] = ...
          append_matrix_memory(p,filesep,INTLAB_MATRIX_MEMORY_USER,INTLAB_MATRIX_FILES_USER);
        setappdata(0,'INTLAB_MATRIX_MEMORY_USER',INTLAB_MATRIX_MEMORY_USER);
        setappdata(0,'INTLAB_MATRIX_FILES_USER',INTLAB_MATRIX_FILES_USER);
        
      end

      % check whether path p belongs to INTLAB path
      if findstr(p,INTLABPath)==1
        
        disp(['adding path ' p ' to matrix memory for INTLAB path'])
        [INTLAB_MATRIX_MEMORY,INTLAB_MATRIX_FILES] = ...
          append_matrix_memory(p,filesep,INTLAB_MATRIX_MEMORY,INTLAB_MATRIX_FILES);
        setappdata(0,'INTLAB_MATRIX_MEMORY',INTLAB_MATRIX_MEMORY);   
        setappdata(0,'INTLAB_MATRIX_FILES',INTLAB_MATRIX_FILES);

        % check whether @-directory exists
        j = findstr(p,filesep);
        p_ = [p filesep '@' p(j(end)+1:end)];
        if exist(p_,'dir')
          disp(['adding path ' p_ ' to matrix memory for complete path'])
          [INTLAB_MATRIX_MEMORY,INTLAB_MATRIX_FILES] = ...
            append_matrix_memory(p_,filesep,INTLAB_MATRIX_MEMORY,INTLAB_MATRIX_FILES);
          setappdata(0,'INTLAB_MATRIX_MEMORY',INTLAB_MATRIX_MEMORY);
          setappdata(0,'INTLAB_MATRIX_FILES',INTLAB_MATRIX_FILES);
        end
        
      end
            
      P = P(i+1:end);
      
    end
    
    cd(dirold)
    disp([int2str(size(INTLAB_MATRIX_FILES,1)) ' different names found in INTLAB paths'])
    disp([int2str(size(INTLAB_MATRIX_FILES_USER,1)) ' different names found in user-defined paths'])
    disp([int2str(size(INTLAB_MATRIX_FILES_ALL,1)) ' different names found in all paths'])
    save(intvalinit('helppfile'), ...
       'INTLAB_MATRIX_MEMORY','INTLAB_MATRIX_FILES', ...
       'INTLAB_MATRIX_MEMORY_USER','INTLAB_MATRIX_FILES_USER', ...
       'INTLAB_MATRIX_MEMORY_ALL','INTLAB_MATRIX_FILES_ALL', ...
       '-mat')
    return
    
  end
  
  % find name
  file = intvalinit('helppfile');
  if ~exist(file,'file')
    error('helpp file not yet generated. Call "helpp" without argument first.')
  end
  load(file)
  defaultlevel = 6;
  minlevel = 3;
  delta = 3;                            % if level not specified, levels from maxlevel-delta on are displayed
  name = lower(name);
  levelspecified = ( nargin==3 );
  
  if nargin==1                                  % input only file name, no extra parameter
    MEMORY = getappdata(0,'INTLAB_MATRIX_MEMORY');
    FILES = getappdata(0,'INTLAB_MATRIX_FILES');
    level = defaultlevel - max(0,5-length(name));  % level smaller for short names
  else                                          % at least two parameters
    all = lower(all);
    if isequal(all,'all')
      MEMORY = getappdata(0,'INTLAB_MATRIX_MEMORY_ALL); 
      FILES = getappdata(0,'INTLAB_MATRIX_FILES_ALL);
      if nargin~=3
        level = defaultlevel - max(0,5-length(name));  % level smaller for short names
      end
    elseif isequal(all,'user')
      MEMORY = getappdata(0,'INTLAB_MATRIX_MEMORY_USER); 
      FILES = getappdata(0,'INTLAB_MATRIX_FILES_USER);
      if nargin~=3
        level = defaultlevel - max(0,6-length(name));  % level smaller for short names
      end
    else
      MEMORY = getappdata(0,'INTLAB_MATRIX_MEMORY); 
      FILES = getappdata(0,'INTLAB_MATRIX_FILES);
      level = all;
      levelspecified = 1;
    end
  end
  
  maxlen = size(MEMORY,2);
  if length(name)>maxlen                        % cut name to maxlen    
    name = name(1:maxlen);
  end
  
  if isempty(findstr(name,'.'))                 % add file extension .m if not specified
    name = [name '.m'];
  end
  
  name = [deblank(name) ' '];                   % exactly one trailing blank
  d_name = double(name);                        % convert to double
  k = 256*(d_name(1:end-1))+d_name(2:end);      % index of name
  Level = sum(MEMORY(:,k),2);
  Index = find( Level >= level );               % minimum concordance
  if any(Index(:)) & ( ~levelspecified )        % display only levels not too far from maximum
    maxlevel = max(Level(Index));
    Index = find( Level >= max(maxlevel-delta,minlevel));
  end
  ll = sum(MEMORY(Index,k),2)/length(name);
  [sortLevel,I] = sort(Level(Index));
  len = length(sortLevel);
  Names = FILES(Index,:);                       % selected rows corresponding to name
  [ repmat('level ',len,1) int2str(full(sortLevel)) repmat('  ',len,1) Names(I,:)]

  
  function [MEMORY,FILES] = append_matrix_memory(path,filesep,MEMORY,FILES)
    cd(path);                       % change path
    D = dir;                        % build matrix memory
    D = D(3:end);
    files = [];
    if ~isempty(D)
      j = 0;
      for i=1:length(D)
        if ~D(i).isdir
          j = j+1;
          files{j} = D(i).name;                     % individual file name in directory path       
        end
      end
      if ~isempty(files)
        names = double(strvcat(lower(files)));      % lower case filenames converted to double for memory index
        Names = double(strvcat(files));             % filenames converted to double for later display
        [m,n] = size(Names);                        % m  number of filenames; n  max length of filename
        Index = 256*names(:,1:end-1) + names(:,2:end);      % array m x (n-1)
        [M,N] = size(MEMORY);
        MEMORY(M+1:M+m,:) = sparse(repmat((1:m)',n-1,1),Index(:),1,m,65536);   % matrix memory 
        FILES = strvcat(FILES,[repmat(path,m,1) repmat(filesep,m,1) Names]);
      end
    end  
