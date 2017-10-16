processingMode  = 1;  % 0: compression, 1: decompression, 2: JP2/KLB format conversion
compressionType = 0;  % 0: KLB file format, 1: JP2 file format

deleteOld       = 1;  % set to "1" to delete source files after processing
maxIterations   = 50; % number of attempts at reading/writing a file in the event of an I/O error before skipping that file

poolWorkers     = 12; % use "0" to enable automated detection of available CPU cores
                      % use "-1" to disable parallel proecessing

%% main loop

switch compressionType
    case 0
        compressionExtension = 'klb';
    case 1
        compressionExtension = 'jp2';
end;

if processingMode == 2
    compressionExtension1 = 'jp2';
    compressionExtension2 = 'klb';
end;

directoryList = cell(1, 1);
directoryList{1} = cd;

if poolWorkers == 0
    poolWorkers = feature('numcores');
    disp(' ');
    disp([num2str(poolWorkers) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

if poolWorkers == -1
    while ~isempty(directoryList)
        currentDir = directoryList{1};
        disp(['processing directory ' currentDir]);
        
        currentSubs = dir(currentDir);
        for i = 1:length(currentSubs)
            if currentSubs(i).isdir && ~strcmp(currentSubs(i).name, '.') && ~strcmp(currentSubs(i).name, '..')
                directoryList = cat(1, directoryList, {[currentDir '\' currentSubs(i).name]});
            end;
        end;
        
        switch processingMode
            case 0
                currentUncompressedFiles = dir([currentDir '\*.tif']);
                
                if ~isempty(currentUncompressedFiles)
                    for currentIndex = 1:numel(currentUncompressedFiles)
                        inputName   = [currentDir '\' currentUncompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentUncompressedFiles(currentIndex).name(1:(end - 3)) compressionExtension];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
                
            case 1
                currentCompressedFiles = dir([currentDir '\*.' compressionExtension]);
                
                if ~isempty(currentCompressedFiles)
                    for currentIndex = 1:numel(currentCompressedFiles)
                        inputName   = [currentDir '\' currentCompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentCompressedFiles(currentIndex).name(1:(end - 3)) 'tif'];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
                
            case 2
                currentCompressedFiles = dir([currentDir '\*.' compressionExtension1]);
                
                if ~isempty(currentCompressedFiles)
                    for currentIndex = 1:numel(currentCompressedFiles)
                        inputName   = [currentDir '\' currentCompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentCompressedFiles(currentIndex).name(1:(end - 3)) compressionExtension2];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
        end;
        
        directoryList(1) = [];
    end;
else
    disp(' ');
    if matlabpool('size') > 0
        matlabpool('close');
    end;
    matlabpool(poolWorkers);
    disp(' ');
    
    warning off;
    
    while ~isempty(directoryList)
        currentDir = directoryList{1};
        disp(['processing directory ' currentDir]);
        
        currentSubs = dir(currentDir);
        for i = 1:length(currentSubs)
            if currentSubs(i).isdir && ~strcmp(currentSubs(i).name, '.') && ~strcmp(currentSubs(i).name, '..')
                directoryList = cat(1, directoryList, {[currentDir '\' currentSubs(i).name]});
            end;
        end;
        
        switch processingMode
            case 0
                currentUncompressedFiles = dir([currentDir '\*.tif']);
                
                if ~isempty(currentUncompressedFiles)
                    parfor currentIndex = 1:numel(currentUncompressedFiles)
                        inputName   = [currentDir '\' currentUncompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentUncompressedFiles(currentIndex).name(1:(end - 3)) compressionExtension];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
                
            case 1
                currentCompressedFiles = dir([currentDir '\*.' compressionExtension]);
                
                if ~isempty(currentCompressedFiles)
                    parfor currentIndex = 1:numel(currentCompressedFiles)
                        inputName   = [currentDir '\' currentCompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentCompressedFiles(currentIndex).name(1:(end - 3)) 'tif'];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
                
            case 2
                currentCompressedFiles = dir([currentDir '\*.' compressionExtension1]);
                
                if ~isempty(currentCompressedFiles)
                    parfor currentIndex = 1:numel(currentCompressedFiles)
                        inputName   = [currentDir '\' currentCompressedFiles(currentIndex).name];
                        outputName  = [currentDir '\' currentCompressedFiles(currentIndex).name(1:(end - 3)) compressionExtension2];
                        fileRead    = 0;
                        fileWritten = 0;
                        
                        iterations = 1;
                        while ~fileRead && (iterations <= maxIterations)
                            try
                                stack = readImage(inputName);
                                fileRead = 1;
                            catch errorMessage
                                disp(['error reading file ' inputName]);
                                pause(1);
                                iterations = iterations + 1;
                            end;
                        end;
                        if fileRead
                            iterations = 1;
                            while ~fileWritten && (iterations <= maxIterations)
                                try
                                    writeImage(stack, outputName);
                                    fileWritten = 1;
                                catch errorMessage
                                    disp(['error writing file ' outputName]);
                                    pause(1);
                                    iterations = iterations + 1;
                                end;
                            end;
                        end;
                        if ~fileWritten
                            disp(['skipping file ' inputName]);
                        end;
                        if deleteOld && fileWritten
                            delete(inputName);
                        end;
                    end;
                end;
        end;
        
        directoryList(1) = [];
    end;
    
    warning on;
    
    disp(' ');
    if matlabpool('size') > 0
        matlabpool('close');
    end;
    disp(' ');
end;