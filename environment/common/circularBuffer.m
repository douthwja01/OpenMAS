%% CIRCULAR BUFFER (circularBuffer.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class provides the functionality of a circular buffer to provide
% sequential overriding of allocated memory.

classdef circularBuffer    
    properties
        data;       % Container for the numeric values
    end
    % Methods on the circular buffer
    methods
        % Constructor procedure
        function obj = circularBuffer(data)
            obj.data = data; % Just store, regardless of its type
        end
        % Reference procedure
        function out = subsref(obj,s)
            
            switch s(1).type
                case '()'
                    sizePerDim = size(obj.data);                           % Get the buffer dimensions
                    dimSelect = repmat({':'},1,numel(sizePerDim));         % Default to all on unselected dimensions
                    % Accessed as a vector
                    if numel(s(1).subs) == 1
                        % Accessed as a vector
                        [~,b] = max(sizePerDim);
                        dimSelect{b} = s(1).subs{1};
                        refDim = dimSelect;
                    else
                        % Accessed as a matrix
                        refDim = s(1).subs;
                    end
                    % Move through each input dimension and assign the
                    % input to that dimension
                    newDim = cell(size(refDim));
                    for i = 1:numel(refDim)                                % For each of the requested dimensions
                        if ~strcmp([refDim{i}],':')                        
                            newDim{i} = mod(refDim{i}-1,sizePerDim(i)) + 1;% if the dimension is numeric
                        else
                            newDim{i} = refDim{i};
                        end
                    end
                    out = obj.data(newDim{:});
                case '.'
                    % Behave as if the "." operator is standard
                    out = builtin('subsref',obj,s);
                otherwise
                    error(['Error(circularBuffer): Use bracket always... Currently used ',s(1).type,' !!!']);
            end
        end
        % Assignment procedure
        function obj = subsasgn(obj,s,input)
                 
            % where ":"
            % inputSize == bufferSize
            % where "2"
            
            % Input sanity check
%             inputSize = size(input);
%             bufferSize = size(obj.data);
%             for i = 1:numel(s(1).subs)
%                 % Where ":" is assigned, dimensions must match
%                 if strcmp([s(1).subs{i}],':')
%                     Astr = sprintf('Dimension %d: Buffer dimension is %d, where Input dimension is %d.',i,bufferSize(i),inputSize(i));
%                     assert(inputSize(i) == bufferSize(i),Astr);
%                 end
%             end

            % Type of sub-assign
            switch s(1).type
                case '()'                   
                    % Get the buffer dimensions
                    sizePerDim = size(obj.data);
                    assignmentDim = s(1).subs;
                    
                    % Default behaviour when incorrect dimensions are given
%                     dimSelect = repmat({':'},1,numel(sizePerDim));
                    dimSelect = repmat({1},1,numel(sizePerDim));           
                    % We cannot allow allocation of one value to entire row
                    
                    % Accessed as a vector
                    if numel(assignmentDim) == 1 
                        % Accessed as a vector
                        [~,b] = max(sizePerDim);
                        dimSelect{b} = assignmentDim{1};
                        assignmentDim = dimSelect;
                    else
                        % Accessed as a matrix
                        assignmentDim = s(1).subs;
                    end
                                      
                    % if accessed as var(5):
                    % - Define var(:,5) if row
                    % - Define var(5,:) if column
                    
                    % Move through each input dimension and assign the
                    % input to that dimension
                    for i = 1:numel(assignmentDim)
                        % For each of the requested dimensions
                        if ~strcmp([assignmentDim{i}],':')
                            % if the dimension is numeric
                            newDim{i} = mod(assignmentDim{i}-1,sizePerDim(i)) + 1;
                        else
                            % All entries
                            newDim{i} = assignmentDim{i};
                        end
                    end
                    % Attempt to assign new value
                    try
                        obj.data(newDim{:}) = input;                           % Assign new value to buffer
                    catch insertError
                        warning('Unable to insert value into circular buffer, do the dimensions match?');
                        rethrow(insertError);
                    end
                otherwise
                    error(['Error(circularBuffer): Use bracket always... Currently used ',s(1).type,' !!!']);
            end
        end
        
        % Unpack procedure
        function [dataSequence] = unpack(obj,unpackDim,ind0)
            % This function will unpack a provided buffer along the
            % dimension provided given an initial starting position
            % startInd.
            % INPUTS:
            % unpackDim - The dimension along which the data is resequenced
            % ind0      - The reference (or position of last entry).
            % OUTPUT:
            % dataSeq   - The buffer data sequenced with ind0 first
            
            % Dimensions of the buffer
            bufferDim = size(obj.data);
            
            % Assume we are only moving along the second axis
            assert(numel(unpackDim) == 1 && unpackDim <= numel(bufferDim),'The unpack dimension is invalid.');
            assert(isnumeric(unpackDim) && isnumeric(ind0),'Dimension specifier must be numeric.');
            
            ind0 = double(ind0); % To remove issues with unsigned indicies
            
            % Default dimension selection
            dimSelectA = repmat({':'},1,numel(bufferDim));
            dimSelectB = dimSelectA;
            dataSequence = obj.data;        % Container of same type and size
            
            % while t is less than the length of the buffer in the unpack
            % dimension
            t = 1;                          % Steps back
            while t < (bufferDim(unpackDim) + 1)
                bInd = mod((ind0 - t),bufferDim(unpackDim)) + 1;
                % Modify data based on the selected unpacking dimension
                dimSelectA{unpackDim} = t;
                dimSelectB{unpackDim} = bInd;
                dataSequence(dimSelectA{:}) = obj.data(dimSelectB{:,:});
                t = t + 1;                  % Iterate index of new array
            end
        end
    end
end
    
%     methods (Static)
%         % Override entries in a buffer
%         function [data] = export(buffer,ind)
%             % Get requested data from the buffer
%             bufferDim = size(buffer);
%             inputDim  = numel(ind);
%             % Define buffer-relative indicies
%             newDim = {};
%             for i = 1:inputDim
%                 if strcmp(inputDim{i})
%                     % All entries of that dimension
%                     newDim{i} = ':';
%                 else
%                     % Get the equivalent circular index
%                     newDim{i} = mod(inputDim{i}-1,bufferDim(i)) + 1;
%                 end
%             end
%             % Select the data from the new indicies vector
%             data = buffer(newDim{:});
%         end
%         % Extract entries from a buffer
%         function [data] = import(bufferData,ind,input)
%             % Get requested data from the buffer
%             bufferDim = size(bufferData);
%             inputDim  = numel(ind);
%             % Define buffer-relative indicies
%             newDim = {};
%             for i = 1:inputDim
%                 if strcmp(inputDim{i})
%                     % All entries of that dimension
%                     newDim{i} = ':';
%                 else
%                     % Get the equivalent circular index
%                     newDim{i} = mod(inputDim{i}-1,bufferDim(i)) + 1;
%                 end
%             end
%             % Override the buffer entry with the input
%             bufferData(newDim{:}) = input;
%             data = bufferData;
%         end