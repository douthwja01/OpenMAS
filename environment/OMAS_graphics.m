%% OMAS GRAPHICS (OMAS_graphics.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This utility provides tools for interacting with 3D models within OMAS
% based on the 'stlTools' library. This function does in no way take credit
% for the design of the functions inclosed, only for their arrangment as
% part of the OMAS package.

% Author: James A. Douthwaite 08/03/18

classdef OMAS_graphics
    %% GENERAL GEOMETRY MANIPULATION 
    methods (Static, Access = public)
        % SORT POINT CLOUD IN A CLOCK WISE ABOUT A VECTOR  
        function [verticesCW] = sortVerticesCW(planarNormal,centroid,vertexData)
            % This function takes a list of vertex data and sorts them into
            % angularly in relation to the planar normal vector.
            
            % INPUT HANDLING
            if numel(centroid) == 2
                centroid = [centroid;0];
            end
            % THE REFERENCE
            unit_reference = [1;0;0];

            % DEFINE THE VECTORS RELATIVE TO THE CENTROID
            vertexVectors = vertexData - centroid';
            
            tempMatrix = zeros(size(vertexVectors,1),2);
            for i = 1:size(vertexVectors,1)
                unit_v = vertexVectors(i,:)/norm(vertexVectors(i,:));
                vertexNormal = cross(unit_reference,unit_v);
                % THE ROTATION DIRECTION
                tempMatrix(i,1) = -sign(dot(vertexNormal,planarNormal));
                % THE ROTATION MAGNITUDE
                tempMatrix(i,2) = acos(dot(unit_reference,unit_v));
                % MAP TO CCW ROTATION
                if tempMatrix(i,1) < 0
                    tempMatrix(i,2) = 2*pi - tempMatrix(i,2) ;
                end
            end
            % SORT THE VERTICE BY ANGLE RELATIVE 
            [~,sortedIndices] = sortrows(tempMatrix(:,2),1);
            % GET THE SORTED VERTEX ARRAY
            verticesCW = vertexData(sortedIndices,:);
        end
        % SORT POINT CLOUD IN A COUNTER-CLOCK WISE ABOUT A VECTOR  
        function [verticesCCW] = sortVerticesCCW(planarNormal,centroid,vertexData)
            % This function takes a list of vertex data and sorts them into
            % angularly in relation to the planar normal vector.
            
            % INPUT HANDLING
            if numel(centroid) < 3
                centroid = [centroid;0];
            end
            % THE REFERENCE
            unit_reference = [1;0;0];

            % DEFINE THE VECTORS RELATIVE TO THE CENTROID
            vertexVectors = vertexData - centroid';
            
            tempMatrix = zeros(size(vertexVectors,1),2);
            for i = 1:size(vertexVectors,1)
                unit_v = vertexVectors(i,:)/norm(vertexVectors(i,:));
                vertexNormal = cross(unit_reference,unit_v);
                % THE ROTATION DIRECTION
                tempMatrix(i,1) = sign(dot(vertexNormal,planarNormal));
                % THE ROTATION MAGNITUDE
                tempMatrix(i,2) = acos(dot(unit_reference,unit_v));
                % MAP TO CCW ROTATION
                if tempMatrix(i,1) < 0
                    tempMatrix(i,2) = 2*pi - tempMatrix(i,2) ;
                end
            end
            % SORT THE VERTICE BY ANGLE RELATIVE 
            [~,sortedIndices] = sortrows(tempMatrix(:,2),1);
            % GET THE SORTED VERTEX ARRAY
            verticesCCW = vertexData(sortedIndices,:);
            
%             % DEBUG PLOTS
%             figure(3)
%             hold on; axis equal; grid on;
%             scale = 10;
%             q = quiver3(centroid(1),centroid(2),centroid(3),...
%                         scale*unit_reference(1),scale*unit_reference(2),scale*unit_reference(3),'r');
%             q = quiver3(centroid(1),centroid(2),centroid(3),...
%                         scale*planarNormal(1),scale*planarNormal(2),scale*planarNormal(3),'g');
%             for i = 1:size(verticesCCW,1)
%                 plot3(verticesCCW(i,1),verticesCCW(i,2),verticesCCW(i,3),'ro');
%             end
        end
        % PLANAR PROHECTION OF A 3D GEOMETRY
        function [planarVertices] = geometryPlanarProjection(planarNormal,centroid,vertexData)
            % This function returns a vertex projections on a plane defined
            % by the normal vector provided
            
            % INPUT HANDLING
            assert(size(planarNormal,1) == size(centroid,1),'The planar normal and centroid must be [3x1]');
            assert(size(planarNormal,2) == 1,'The planar normal and centroid must be [3x1]');
            assert(size(vertexData,2) == 3,'The geometry data must be [:x3].');
                             
            % DEFINE THE VECTORS BETWEEN THE CENTROID AND VERTICES
            vertexData = vertexData - centroid';
            planarVertices = zeros(size(vertexData,1),3);
            for i = 1:size(vertexData,1) 
                % The vertex planar projections
                [projection] = OMAS_geometry.vectorPlanarProjection(...
                                    planarNormal,...
                                    vertexData(i,:)');
                % REFORMAT THE PROJECTIONS
                planarVertices(i,:) = projection' + centroid';
            end
            % REMOVE DUPLICATES
            [planarVertices,~] =  unique(planarVertices, 'rows');
            % We cannot gaurantee the geometry will be the same one the
            % unique command is ran. The face data may reference vertices
            % that have been moved or deleted.
        end
        % DEFINE SCALE OF GEOMETRY STRUCTURE
        function [geometry] = scale(geometry,scale)
            % This function scales a geometry in accordance to a provided
            % scalar or vector of dimensional scalars. A geometry is
            % defined as a vertices, faces structure.
            
            assert(numel(scale) == 1 || numel(scale) == 3,'Please provide a valid scale value.');
            % Apply the scaling either by dimension or unilaterally.
            if numel(scale) == 3
                geometry.vertices = geometry.vertices*diag(scale);
            else
                geometry.vertices = geometry.vertices*(eye(3)*scale);
            end
        end
        % NORMALISE THE STL FILE
        function [geometry] = normalise(geometry)
            % This function normalises the STL geometry to allow it to be
            % scaled appropriately.
            % GEOMETRY CHECK
            assert(size(geometry.vertices,2) == 3 && size(geometry.faces,2) == 3,...
                'The provided patch object has invalid vertex and faces assignments.');
            
            % GET THE VERTEX MAGNITUDES RELATIVE TO THE ORIGIN
            verticesNorms = abs(geometry.vertices);
            % MAXIMUMS IN THE FIRST AND SECOND DIMENSIONS
            dimMaximal = max(max(abs(verticesNorms),[],1),[],2);
            % SCALE GEOMETRY SPECIFIED BY VIRTUAL.radius
            geometry.vertices = (geometry.vertices/dimMaximal);
        end
        % REMOVE FACE/VERTEX DUPLICATES
        function [geometry] = removeDuplicateVertices(geometry)
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
            
            if ~isfield(geometry,'vertices')
                error('The geometry does not have a specified "vertex" list.');
            end
            if ~isfield(geometry,'faces')
                error('The geometry does not have a specified "face" triangulation list.');
            end
            
            % REMOVE DUPLICATES
            [geometry.vertices,~,indexn] =  unique(geometry.vertices, 'rows');
            geometry.faces = indexn(geometry.faces);
        end
        % CALCULATE GEOMETRY SURFACE NORMALS
        function [normals]  = normals(geometry)
            % This function computes the surface normals of a defined
            % geometry structure
            normals = zeros(size(geometry.faces));
            for face = 1:size(geometry.faces,1)
                % MEMBERS IF THE PLANE
                memberID_A = geometry.faces(face,1);
                memberID_B = geometry.faces(face,2);
                memberID_C = geometry.faces(face,3);
                % SURFACE DEFINING VECTORS
                BA = geometry.vertices(memberID_B,:) - geometry.vertices(memberID_A,:);
                CA = geometry.vertices(memberID_C,:) - geometry.vertices(memberID_A,:);
                % THE NORMAL
                normals(face,:) = cross(BA,CA);
                normals(face,:) = normals(face,:)/norm(normals(face,:));
            end
        end
    end
    
    %% UNIVERSAL DRAWING MECHANISMS
    methods (Static)
        % GET HIT-BOX GEOMETRY
        function [hitBoxGeometry] = getHitBoxGeometry(VIRTUAL,geometry)
            % Input check
            assert(isa(VIRTUAL.hitBoxType,'uint8'),'Object hit-box type must be of type "uint8".');
            % Derive hit-box geometry
            switch VIRTUAL.hitBoxType
                case OMAS_hitBoxType.none
                    % Assemble the spherical constraint volume
                    hitBoxGeometry = [];
                case OMAS_hitBoxType.spherical
                    % Assemble the spherical constraint volume
                    hitBoxGeometry = OMAS_graphics.defineSphere(zeros(3,1),VIRTUAL.radius,10);
                case OMAS_hitBoxType.AABB
                    % Assemble the AABB constraint volume
                    R = OMAS_geometry.quaternionToRotationMatrix(VIRTUAL.quaternion); 
                    % Rotate the body before defining constraint volume
                    geometry.vertices = R*geometry.vertices;
                    % DEFINE AN ALIGNED AABB CUBOID
                    hitBoxGeometry = OMAS_graphics.defineCuboid(min(geometry.vertices),...
                                                                  max(geometry.vertices));
                case OMAS_hitBoxType.OBB
                    % DEFINE AN ALIGNED OBB CUBOID
                    hitBoxGeometry = OMAS_graphics.defineCuboid(min(geometry.vertices),...
                                                                  max(geometry.vertices));                    
                    % Assemble the OBB constraint volume
                    R = OMAS_geometry.quaternionToRotationMatrix(VIRTUAL.quaternion); 
                    % Rotate the hit-box to be aligned with the geometry
                    hitBoxGeometry.vertices = R*hitBoxGeometry.vertices;
                otherwise
                    error('Hit-box type not recognised');
            end
        end
        % DEFINE CUBOID (MINOR) FROM RADIUS
        function [geometry,minExtents,maxExtents] = defineCuboidFromRadius(center,radius)
            % This function creates a set of vertices for a cube
            % encapsuated in a sphere of a given radius.
            assert(numel(center) == 3,'The cuboid center must be a cartesian vector.');
            assert(numel(radius) == 1,'The radius of the sphere must be a scalar and non-zero.');
            % RATE OF DEMENSIONAL EXPANSION
            h = radius/1.7321;
            % DEFINE THE CUBOID EXTENTS
            minExtents = center - h;
            maxExtents = center + h;
            % DEFINE THE CUBOID VERTICES
            [geometry] = OMAS_graphics.defineCuboid(minExtents,maxExtents);
        end
        % DEFINE CUBOID
        function [geometry,minExtents,maxExtents] = defineCuboid(minExtents,maxExtents)
            % Return a matrix of point defining a cuboid scaled to that of
            % a dimensions provided.
            
            % Define vertex data from limits
            vertices = zeros(8,3);
            vertices(1,:) = [maxExtents(1),maxExtents(2),maxExtents(3)];
            vertices(2,:) = [maxExtents(1),maxExtents(2),minExtents(3)];
            vertices(3,:) = [maxExtents(1),minExtents(2),minExtents(3)];
            vertices(4,:) = [maxExtents(1),minExtents(2),maxExtents(3)];
            vertices(5,:) = [minExtents(1),minExtents(2),minExtents(3)];
            vertices(6,:) = [minExtents(1),minExtents(2),maxExtents(3)];
            vertices(7,:) = [minExtents(1),maxExtents(2),maxExtents(3)];
            vertices(8,:) = [minExtents(1),maxExtents(2),minExtents(3)];
            geometry.vertices = vertices;
            % Define face connectivity matrix
            geometry.faces =  [  1     2     7
                                 1     4     2
                                 1     7     4
                                 2     3     8
                                 2     4     3
                                 2     8     7
                                 3     4     6
                                 3     5     8
                                 3     6     5
                                 4     7     6
                                 5     6     8
                                 6     7     8];
        end
        % DRAW SPHERE
        function [geometry] = defineSphere(position,radius,faces)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius. 'figureHandle' is used
            % to provide the function context.
            
            % INPUT CHECKING
            if nargin < 3
                faces = 10;
            end
            % DEFINE SPHERE TRIANGULATION
            [X,Y,Z] = sphere(faces);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
            % CONVERT TO PATCH OBJECT
            [geometry.faces,geometry.vertices,~] = surf2patch(X,Y,Z,'triangles');
            % DEFINE CENTROID
            geometry.centroid = position'; %reducepatch(P, R)
        end
        % DRAW UNIT TRIAD
        function [triadHandle] = drawTriad(figureHandle,position,R,scale)
            % Draw a unit triad at a cartesian position, rotated by R.
            % FigureHandl Here is used to bring the current figure to the
            % function context.
            
            % INPUT CHECL
            if nargin < 4
                scale = 1;
            end
            
            colourVector = 'rgb';
            triadVectors = scale*eye(3);
            for axis = 1:size(triadVectors,2)
                triadVectors(:,axis) = R*triadVectors(:,axis);
                % Draw vectors
                triadHandle(axis) = quiver3(gca,position(1),position(2),position(3),triadVectors(1,axis),triadVectors(2,axis),triadVectors(3,axis),colourVector(axis));
            end
        end
    end

    %% STL IMPORTATION  
    methods (Static, Access = public)
        % PREPARE THE PATCH FILES FOR VISUALS
        function [geometry,successFlag] = importStlFromFile(filename)
            % This function prepares object STL files for presentation. The assumption
            % is that object stl file has the same name as the object being simulated.
            
            try
                % GET THE STL FILE
                geometry = OMAS_graphics.stlRead(filename);
                successFlag = 1;
            catch
                % NO STL WAS FOUND BY THAT FILE NAME (fail quietly)
                geometry = [];
                successFlag = 0;
            end
            % ENSURE UNIQUE VERTICES
            if successFlag
                geometry = OMAS_graphics.removeDuplicateVertices(geometry);
            end
        end
        % READ AND STL INTO VERTICIES, FACES
        function [geometry, name] = stlRead(fileName)
            %STLREAD reads any STL file not depending on its format
            %V are the vertices
            %F are the faces
            %N are the normals
            %NAME is the name of the STL object (NOT the name of the STL file)
            
            [format,isSuccessful] = OMAS_graphics.stlGetFormat(fileName);
            
            if strcmp(format,'ascii')
                [geometry,name] = OMAS_graphics.stlReadAscii(fileName);
            elseif strcmp(format,'binary')
                [geometry,name] = OMAS_graphics.stlReadBinary(fileName);
            end
        end
        % IDENTIFY STL TYPE
        function [format,isSuccessful] = stlGetFormat(fileName)
            %STLGETFORMAT identifies the format of the STL file and returns 'binary' or
            %'ascii'
            
            fid = fopen(fileName);
            if fid == -1        % Unsuccessful file load
                isSuccessful = 0;
            else
                isSuccessful = 1;
            end
            
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
                fseek(fid,0,-1);                                           % go to the beginning of the file
                header = strtrim(char(fread(fid,80,'uchar')'));            % trim leading and trailing spaces
                isSolid = strcmp(header(1:min(5,length(header))),'solid'); % take first 5 char
                fseek(fid,-80,1);                                          % go to the end of the file minus 80 characters
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
        function [geometry, name] = stlReadBinary(fileName)
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
            % Geometric structure
            geometry = struct('vertices',v,'faces',f,'normals',n);
            [geometry] = OMAS_graphics.removeDuplicateVertices(geometry);
        end
        % INTERPRET AN ASCII STL FILE
        function [geometry, n, name] = stlReadAscii(fileName)
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
            cellcontent = textscan(fid,'%s','delimiter','\n');             % read all the file and put content in cells
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
            n = str2double(normals(:,13:end));
            
            % read the vertex coordinates (vertices)
            vertices = char(content(logical(strncmp(content,'vertex',6))));
            v = str2double(vertices(:,7:end));
            nvert = size(vertices,1);                       % number of vertices
            nfaces = sum(strcmp(content,'endfacet'));       % number of faces
            if (nvert == 3*nfaces)
                f = reshape(1:nvert,[3 nfaces])';           % create faces
            end
            % Geometric structure
            geometry = struct('vertices',v,'faces',f,'normals',n);
            % slim the file (delete duplicated vertices)
            [geometry] = OMAS_graphics.removeDuplicateVertices(geometry);
        end
        % PLOT THE STL
        function stlPlot(geometry, name)
            %STLPLOT is an easy way to plot an STL object
            %V is the Nx3 array of vertices
            %F is the Mx3 array of faces
            %NAME is the name of the object, that will be displayed as a title
            
            figure;
            patch(geometry,...
                'FaceColor',[0.8 0.8 1.0], ...
                'EdgeColor','none',...
                'FaceLighting','gouraud',...
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
end