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
        % Reference proceedure
        function out = subsref(obj,s)
            
            switch s(1).type
                case '()'
                    % Get the buffer dimensions
                    bufferDim = size(obj.data);
                    % Get the input dimensions
                    inputDim = s(1).subs;
                    % Move through each input dimension and assign the
                    % input to that dimension
                    for i = 1:numel(inputDim)
                        % For each of the requested dimensions
                        if ~strcmp([inputDim{i}],':')
                            % if the dimension is numeric
                            newDim{i} = mod(inputDim{i}-1,bufferDim(i)) + 1;
                        else
                            % All entries
                            newDim{i} = inputDim{i};
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
            inputSize = size(input);
            bufferSize = size(obj.data);
            for i = 1:numel(s(1).subs)
                % Where ":" is assigned, dimensions must match
                if strcmp([s(1).subs{i}],':')
                    Astr = sprintf('Dimension %d: Buffer dimension is %d, where Input dimension is %d.',i,bufferSize(i),inputSize(i));
                    assert(inputSize(i) == bufferSize(i),Astr);
                end
            end

            % Type of sub-reference
            switch s(1).type
                case '()'
                    % Get the buffer dimensions
                    bufferDim = size(obj.data);
                    % Get the input dimensions
                    inputDim = s(1).subs;
                    
                    % Move through each input dimension and assign the
                    % input to that dimension
                    for i = 1:numel(inputDim)
                        % For each of the requested dimensions
                        if ~strcmp([inputDim{i}],':')
                            % if the dimension is numeric
                            newDim{i} = mod(inputDim{i}-1,bufferDim(i)) + 1;
                        else
                            % All entries
                            newDim{i} = inputDim{i};
                        end
                    end
                    % Convert the double matrix back to a circular buffer
                    obj.data(newDim{:}) = input;
                otherwise
                    error(['Error(circularBuffer): Use bracket always... Currently used ',s(1).type,' !!!']);
            end
        end
        
        % Unpack procedure
        function [dataSeq] = unpack(obj,unpackDim,ind0)
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
            
            dimSelectA = repmat({':'},1,numel(bufferDim));
            dimSelectB = dimSelectA;
            dataSeq = obj.data;     % Container of same type and size
            
            t = 0; % Steps back
            while t < bufferDim(unpackDim)
                tind = ind0 - t;
                bInd = mod(tind - 1,bufferDim(unpackDim)) + 1;
                t = t + 1;
                
                % Modify data based on the selected unpacking dimension
                dimSelectA{unpackDim} = t;
                dimSelectB{unpackDim} = bInd;
                dataSeq(dimSelectA{:}) = obj.data(dimSelectB{:});
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