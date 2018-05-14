%% OMAS VISUAL (OMAS_visualTools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This utility provides tools for interacting with 3D models within OMAS
% based on the 'stlTools' library. This function does in no way take credit
% for the design of the functions inclosed, only for their arrangment as
% part of the OMAS package.

% Author: James A. Douthwaite 08/03/18

classdef OMAS_visualTools
    
    methods (Static, Access = public)
        
        % PREPARE THE PATCH FILES FOR VISUALS
        function [patchObject,successFlag] = getSTLPatchObject(filename)
            % This function prepares object STL files for presentation. The assumption
            % is that object stl file has the same name as the object being simulated.
            
            try
                % GET THE STL FILE
                [v, f, ~, ~] = OMAS_visualTools.stlRead(filename);
                patchObject = struct('vertices',v,'faces',f);
                successFlag = 1;
            catch
                % NO STL WAS FOUND BY THAT FILE NAME (fail quietly)
                patchObject = [];
                successFlag = 0;
            end
        end
        
        % REMOVE FACE/VERTEX DUPLICATES
        function [vnew, fnew]= stlSlimVerts(v, f)
            % PATCHSLIM removes duplicate vertices in surface meshes.
            %
            % This function finds and removes duplicate vertices.
            %
            % USAGE: [v, f]=patchslim(v, f)
            %
            % Where v is the vertex list and f is the face list specifying vertex
            % connectivity.
            %
            % v contains the vertices for all triangles [3*n x 3].
            % f contains the vertex lists defining each triangle face [n x 3].
            %
            % This will reduce the size of typical v matrix by about a factor of 6.
            %
            % For more information see:
            %  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
            %
            % Francis Esmonde-White, May 2010
            
            if ~exist('v','var')
                error('The vertex list (v) must be specified.');
            end
            if ~exist('f','var')
                error('The vertex connectivity of the triangle faces (f) must be specified.');
            end
            
            [vnew, indexm, indexn] =  unique(v, 'rows');
            fnew = indexn(f);
        end
        % PLOT THE STL
        function stlPlot(v, f, name)
            %STLPLOT is an easy way to plot an STL object
            %V is the Nx3 array of vertices
            %F is the Mx3 array of faces
            %NAME is the name of the object, that will be displayed as a title
            
            figure;
            object.vertices = v;
            object.faces = f;
            patch(object,'FaceColor',       [0.8 0.8 1.0], ...
                'EdgeColor',       'none',        ...
                'FaceLighting',    'gouraud',     ...
                'AmbientStrength', 0.15);
            
            % Add a camera light, and tone down the specular highlighting
            camlight('headlight');
            material('dull');
            
            % Fix the axes scaling, and set a nice view angle
            axis('image');
            view([-135 35]);
            grid on;
            title(name);
        end
    end
    %% STL IMPORTATION
    methods (Static, Access = public)
        % READ AND STL INTO VERTICIES, FACES
        function [v, f, n, name] = stlRead(fileName)
            %STLREAD reads any STL file not depending on its format
            %V are the vertices
            %F are the faces
            %N are the normals
            %NAME is the name of the STL object (NOT the name of the STL file)
            
            format = OMAS_visualTools.stlGetFormat(fileName);
            if strcmp(format,'ascii')
                [v,f,n,name] = OMAS_visualTools.stlReadAscii(fileName);
            elseif strcmp(format,'binary')
                [v,f,n,name] = OMAS_visualTools.stlReadBinary(fileName);
            end
        end
        % IDENTIFY STL TYPE
        function format = stlGetFormat(fileName)
            %STLGETFORMAT identifies the format of the STL file and returns 'binary' or
            %'ascii'
            
            fid = fopen(fileName);
%             if fid == -1
%                 warning('Unable to open the file: %s',fileName);
%             end
            
            % Check the file size first, since binary files MUST have a size of 84+(50*n)
            fseek(fid,0,1);         % Go to the end of the file
            fidSIZE = ftell(fid);   % Check the size of the file
            if rem(fidSIZE-84,50) > 0
                format = 'ascii';
            else
                % Files with a size of 84+(50*n), might be either ascii or binary...
                % Read first 80 characters of the file.
                % For an ASCII file, the data should begin immediately (give or take a few
                % blank lines or spaces) and the first word must be 'solid'.
                % For a binary file, the first 80 characters contains the header.
                % It is bad practice to begin the header of a binary file with the word
                % 'solid', so it can be used to identify whether the file is ASCII or
                % binary.
                fseek(fid,0,-1); % go to the beginning of the file
                header = strtrim(char(fread(fid,80,'uchar')')); % trim leading and trailing spaces
                isSolid = strcmp(header(1:min(5,length(header))),'solid'); % take first 5 char
                fseek(fid,-80,1); % go to the end of the file minus 80 characters
                tail = char(fread(fid,80,'uchar')');
                isEndSolid = findstr(tail,'endsolid');
                
                % Double check by reading the last 80 characters of the file.
                % For an ASCII file, the data should end (give or take a few
                % blank lines or spaces) with 'endsolid <object_name>'.
                % If the last 80 characters contains the word 'endsolid' then this
                % confirms that the file is indeed ASCII.
                if isSolid && isEndSolid
                    format = 'ascii';
                else
                    format = 'binary';
                end
            end
            fclose(fid);
        end
        % INTERPRET A BINARY STL FILE
        function [v, f, n, name] = stlReadBinary(fileName)
            %STLREADBINARY reads a STL file written in BINARY format
            %V are the vertices
            %F are the faces
            %N are the normals
            %NAME is the name of the STL object (NOT the name of the STL file)
            
            %=======================
            % STL binary file format
            %=======================
            % Binary STL files have an 84 byte header followed by 50-byte records, each
            % describing a single facet of the mesh.  Technically each facet could be
            % any 2D shape, but that would screw up the 50-byte-per-facet structure, so
            % in practice only triangular facets are used.  The present code ONLY works
            % for meshes composed of triangular facets.
            %
            % HEADER:
            % 80 bytes:  Header text
            % 4 bytes:   (int) The number of facets in the STL mesh
            %
            % DATA:
            % 4 bytes:  (float) normal x
            % 4 bytes:  (float) normal y
            % 4 bytes:  (float) normal z
            % 4 bytes:  (float) vertex1 x
            % 4 bytes:  (float) vertex1 y
            % 4 bytes:  (float) vertex1 z
            % 4 bytes:  (float) vertex2 x
            % 4 bytes:  (float) vertex2 y
            % 4 bytes:  (float) vertex2 z
            % 4 bytes:  (float) vertex3 x
            % 4 bytes:  (float) vertex3 y
            % 4 bytes:  (float) vertex3 z
            % 2 bytes:  Padding to make the data for each facet 50-bytes in length
            %   ...and repeat for next facet...
            
            fid = fopen(fileName);
            header = fread(fid,80,'int8'); % reading header's 80 bytes
            name = deblank(native2unicode(header,'ascii')');
            if isempty(name)
                name = 'Unnamed Object'; % no object name in binary files!
            end
            nfaces = fread(fid,1,'int32');  % reading the number of facets in the stl file (next 4 byters)
            nvert = 3*nfaces; % number of vertices
            % reserve memory for vectors (increase the processing speed)
            n = zeros(nfaces,3);
            v = zeros(nvert,3);
            f = zeros(nfaces,3);
            for i = 1 : nfaces % read the data for each facet
                tmp = fread(fid,3*4,'float'); % read coordinates
                n(i,:) = tmp(1:3); % x,y,z components of the facet's normal vector
                v(3*i-2,:) = tmp(4:6); % x,y,z coordinates of vertex 1
                v(3*i-1,:) = tmp(7:9); % x,y,z coordinates of vertex 2
                v(3*i,:) = tmp(10:12); % x,y,z coordinates of vertex 3
                f(i,:) = [3*i-2 3*i-1 3*i]; % face
                fread(fid,1,'int16'); % Move to the start of the next facet (2 bytes of padding)
            end
            fclose(fid);
            % slim the file (delete duplicated vertices)
            [v,f] = OMAS_visualTools.stlSlimVerts(v,f);
        end
        % INTERPRET AN ASCII STL FILE
        function [v, f, n, name] = stlReadAscii(fileName)
            %STLREADASCII reads a STL file written in ASCII format
            %V are the vertices
            %F are the faces
            %N are the normals
            %NAME is the name of the STL object (NOT the name of the STL file)
            
            %======================
            % STL ascii file format
            %======================
            % ASCII STL files have the following structure.  Technically each facet
            % could be any 2D shape, but in practice only triangular facets tend to be
            % used.  The present code ONLY works for meshes composed of triangular
            % facets.
            %
            % solid object_name
            % facet normal x y z
            %   outer loop
            %     vertex x y z
            %     vertex x y z
            %     vertex x y z
            %   endloop
            % endfacet
            %
            % <Repeat for all facets...>
            %
            % endsolid object_name
            
            fid = fopen(fileName);
            cellcontent = textscan(fid,'%s','delimiter','\n'); % read all the file and put content in cells
            content = cellcontent{:}(logical(~strcmp(cellcontent{:},''))); % remove all blank lines
            fclose(fid);
            
            % read the STL name
            line1 = char(content(1));
            if (size(line1,2) >= 7)
                name = line1(7:end);
            else
                name = 'Unnamed Object';
            end
            
            % read the vector normals
            normals = char(content(logical(strncmp(content,'facet normal',12))));
            n = str2num(normals(:,13:end));
            
            % read the vertex coordinates (vertices)
            vertices = char(content(logical(strncmp(content,'vertex',6))));
            v = str2num(vertices(:,7:end));
            nvert = size(vertices,1); % number of vertices
            nfaces = sum(strcmp(content,'endfacet')); % number of faces
            if (nvert == 3*nfaces)
                f = reshape(1:nvert,[3 nfaces])'; % create faces
            end
            
            % slim the file (delete duplicated vertices)
            [v,f] = OMAS_visualTools.stlSlimVerts(v,f);
        end
    end
end