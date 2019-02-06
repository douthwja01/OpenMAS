
classdef OMAS_modelViewer
    properties
        figureProperties;
        targetObject;
        triadPosition = [-2;-2;-2];  % Position of the model
        triadScale = 0.5;
        R        = eye(3);      % The rotation of the body and triad together
        R_offset = eye(3);      % The offset between the body axis triad and geometry
        %triad_b;
        %triad_g;
    end
    methods
        % Constructor
        function obj = OMAS_modelViewer(varargin)
            % Input sanity check
            if numel(varargin) ~= 1
                error('Please provide a target object class to view');
            end
            % Assume input is an object label
            if ischar(varargin{1})
                obj.targetObject = eval(varargin{1});                      % Instantiate the object
            elseif ismember('objectDefinition',superclasses(varargin{1})) % Check for the object root class
                obj.targetObject = varargin{1};
            else
                error('The variable provided is not a valid OpenMAS object.');
            end
            % Get the default figure properties
            [obj.figureProperties] = obj.getFigureProperties();
        end
        % Set the model position
        function [obj] = setTriadPosition(obj,p)
            % Input sanity check
            assert(numel(p) == 3 && size(p,1) == 3,'Position must be a column vector [3x1]');
            % Set the triad's relative position
            obj.triadPosition = p;
        end
        % Set the model rotation
        function [obj] = setRotation(obj,eulers)
            % Input sanity check
            assert(numel(eulers) == 3,'Euler rotations must be a column vector [3x1]');
            % Get the rotation matrix
            obj.R = OMAS_geometry.eulersToRotationMatrix(eulers);
        end
        % Show scene
        function [figureHandle] = show(obj,objectClass)
            
            % Input checking
            if nargin > 1
                target = objectClass;
            else
                target = obj.targetObject;
            end
            % For clarity
            properties = obj.figureProperties;
            
            % Define plot
            figureHandle = figure('Name',target.name);
            ax = axes(figureHandle);
            hold on;

            % Other plot attributes
            title(target.name,'fontweight',properties.fontWeight,'fontsize',properties.titleFontSize);
            xlabel('x(m)','fontweight',properties.fontWeight,'fontSize',properties.axisFontSize);
            ylabel('y(m)','fontweight',properties.fontWeight,'fontSize',properties.axisFontSize);
            zlabel('z(m)','fontweight',properties.fontWeight,'fontSize',properties.axisFontSize);
            set(ax,'FontSize',properties.axisFontSize,'fontWeight',properties.fontWeight);
            set(ax,'Color',properties.axesColor);
            axis('tight','equal');
            grid on;  box on;
            % AXES-ON LOGICAL
            set(ax,'visible',properties.axesVisible)
            % SET ORIENTATION
            view(properties.orientation);
            % REPRESENT GEOMETRY AS A PATCH
            target.GEOMETRY = OMAS_graphics.normalise(target.GEOMETRY);
            patch(ax,...
                'Vertices',(target.GEOMETRY.vertices*obj.R_offset)*obj.R',...
                'Faces',target.GEOMETRY.faces,...
                'FaceColor',properties.faceColour,...
                'EdgeColor',properties.edgeColour,...
                'EdgeAlpha',properties.edgeAlpha,...
                'FaceAlpha',properties.faceAlpha,...
                'FaceLighting',properties.faceLighting,...
                'LineWidth',0.1);
            
            % ADD ALIGNED TRIAD
            triadHandle = OMAS_graphics.drawTriad(figureHandle,zeros(3,1),obj.R*obj.R_offset,obj.triadScale);
            set(triadHandle(1,:),'lineWidth',1.2);
            set(triadHandle(1,:),'lineStyle','-');
            % set(triadHandle(1,:),'Color','k');
            % ADD GLOBAL TRIAD
            triadHandle = OMAS_graphics.drawTriad(figureHandle,obj.triadPosition,eye(3),obj.triadScale);
            set(triadHandle(1,:),'lineWidth',1.2);
            set(triadHandle(1,:),'lineStyle','-');
            % set(triadHandle(1,:),'Color','k');
        end
    end
    methods (Static)
        % Get figure properties from 'OMAS_figureProperties'
        function [figureProperties] = getFigureProperties()            
            % Get the figure properties structure common to OMAS
            figureProperties = OMAS_figureProperties();
            % Specific default parameters
            figureProperties.backgroundColor = 'w';
            figureProperties.orientation = [0 90]; %[75 44];
            figureProperties.axesVisible = 'off';
            figureProperties.axesColor = 'none';
            figureProperties.edgeAlpha = 0.02;
            figureProperties.edgeColour = 'k';
            figureProperties.edgeWidth = 0.1;
            figureProperties.faceAlpha = 0.05;
            figureProperties.faceColour = 'b';
        end
    end
end