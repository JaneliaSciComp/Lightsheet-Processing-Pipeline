function writeImage(imageStack, filename, varargin)
[filepath, name, extension] = fileparts(filename);
if ~isempty(varargin) && ~isempty(intersect(varargin{end},{'.jp2','.tif','.tiff','.klb'}))
    extension = varargin{end};
    varargin(end) = [];
else
    % determine file extension to select appropriate write module
    extension={extension};
end
% set all global options to default settings
transposeFlag = 0;

for currentExtension=extension
    currentExtension=currentExtension{:};
    currentFilename =  [filepath filesep name currentExtension];
    if currentExtension=='.jp2'
        
        % set all options to default settings
        nThreads = min(feature('numcores'), 8);
        compression = 0;
        
        % read optional input parameters
        if ~isempty(varargin)
            for c = 1:2:length(varargin)
                switch lower(varargin{c})
                    case {'nthreads'}
                        nThreads = varargin{c + 1};
                    case {'compression'}
                        compression = varargin{c + 1};
                    case {'transpose'}
                        transposeFlag = varargin{c + 1};
                    otherwise
                        error(['Invalid optional argument for JP2 writer: ' varargin{c}]);
                end;
            end;
        end;
        
        if transposeFlag
            if ndims(imageStack) == 2
                imageStack = permute(imageStack, [2 1]);
            elseif ndims(imageStack) == 3
                imageStack = permute(imageStack, [2 1 3]);
            elseif ndims(imageStack) == 4
                imageStack = permute(imageStack, [2 1 3 4]);
            end;
        end;
        
        writeJP2stack(imageStack, currentFilename, nThreads, compression);
    end     
    if ismember(currentExtension,{'.tif', '.tiff'})
        
        % read optional input parameters
        if ~isempty(varargin)
            for c = 1:2:length(varargin)
                switch lower(varargin{c})
                    case {'transpose'}
                        transposeFlag = varargin{c + 1};
                    otherwise
                        error(['Invalid optional argument for TIFF writer: ' varargin{c}]);
                end;
            end;
        end;
        
        if transposeFlag
            if ndims(imageStack) == 2
                imageStack = permute(imageStack, [2 1]);
            elseif ndims(imageStack) == 3
                imageStack = permute(imageStack, [2 1 3]);
            elseif ndims(imageStack) == 4
                imageStack = permute(imageStack, [2 1 3 4]);
            end;
        end;
        
        writeTIFstack(imageStack, currentFilename, 2^31);
    end   
    if currentExtension == '.klb'
        
       % set all options to default settings
        nThreads = feature('numcores');
        pixelSize = [];
        blockSize = [];
        compressionType = [];
        metadata = [];
        
        %read optional input parameters
        if ~isempty(varargin)
            for c = 1:2:length(varargin)
                switch lower(varargin{c})
                    case {'numthreads'}
                        nThreads = varargin{c + 1};
                    case {'blocksize'}
                        blockSize = varargin{c + 1};
                    case {'pixelsize'}
                        pixelSize = varargin{c + 1};
                    case {'compressiontype'}
                        compressionType = varargin{c + 1};
                    case {'metadata'}
                        metadata = varargin{c + 1};
                    case {'transpose'}
                        transposeFlag = varargin{c + 1};
                    otherwise
                        error(['Invalid optional argument for KLB writer: ' varargin{c}]);
                end;
            end;
        end;
        
        if transposeFlag
            if ndims(imageStack) == 2
                imageStack = permute(imageStack, [2 1]);
            elseif ndims(imageStack) == 3
                imageStack = permute(imageStack, [2 1 3]);
            elseif ndims(imageStack) == 4
                imageStack = permute(imageStack, [2 1 3 4]);
            end;
        end;
        
        writeKLBstack(imageStack, currentFilename, nThreads, pixelSize, blockSize, compressionType, metadata);
    end  
    if ~ismember(currentExtension,{'.jp2','.tif','.tiff','.klb'})
        error 'Format not recognized.'
    end
end;