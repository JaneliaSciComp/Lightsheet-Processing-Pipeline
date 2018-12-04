function processTimepoint(inputFolder, outputFolder, projectionFolder, globalMaskFolder, specimen, angle, timepoints, cameras, channels, dimensions, ...
    startsLeft, startsTop, widths, heights, startsFront, depths, inputType, outputType, correctTIFF, rotationFlag, ...
    medianRange, percentile, segmentFlag, flipHFlag, flipVFlag, splitting, kernelSize, kernelSigma, scaling, ...
    references, dependents, thresholds, loggingFlag, verbose, backgroundValues, jobMemory, t, memoryEstimate)

% -----------------------------------------------------------------------------------------------
% | Image data post-processing for simultaneous multi-view light sheet microscopy               |
% |                                                                                             |
% | Correction of insensitive pixels on the sCMOS camera chip,                                  |
% | foreground segmention and lossless image compression                                        |
% |                                                                                             |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2015                          |
% | Email: kellerp@janelia.hhmi.org                                                             |
% |                                                                                             |
% | Utilizes optimization modules and functions by Fernando Amat, HHMI/Janelia Research Campus: |
% | readKLBstack.mexw64                                                                         |
% | writeKLBstack.mexw64                                                                        |
% -----------------------------------------------------------------------------------------------

timepoint = timepoints(t);

version = 1.07;

configuration = cell(37, 1);

configuration{1}  = version;     configuration{2}  = timepoint;  configuration{3}  = inputFolder;  configuration{4}  = outputFolder; configuration{5}  = projectionFolder; configuration{6} = globalMaskFolder;
configuration{7}  = specimen;    configuration{8} = angle; configuration{9}  = timepoints; configuration{10}  = cameras;      configuration{11} = channels;     configuration{12} = dimensions;
configuration{13} = startsLeft;  configuration{14} = startsTop;  configuration{15} = widths;       configuration{16} = heights;
configuration{17} = startsFront; configuration{18} = depths;
configuration{19} = inputType;   configuration{20} = outputType; configuration{21} = correctTIFF;  configuration{22} = rotationFlag; configuration{23} = medianRange;      configuration{24} = percentile;
configuration{25} = segmentFlag; configuration{26} = flipHFlag;  configuration{27} = flipVFlag;    configuration{28} = splitting;    configuration{29} = kernelSize;       configuration{30} = kernelSigma;
configuration{31} = scaling;     configuration{32} = references; configuration{33} = dependents;   configuration{34} = thresholds;
configuration{35} = loggingFlag; configuration{36} = verbose;    configuration{37} = backgroundValues;
configuration{38} = [jobMemory(1) memoryEstimate];

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

if inputType == 4
    if verbose
        disp(['Processing sequence ' num2str(timepoint, '%.4d')]);
    end;
    
    save([outputFolder filesep 'Configuration.mat'], 'configuration');
    
    for h = 1:numel(channels)
        source = [inputFolder filesep 'ch' num2str(channels(h)) '.xml'];
        target = [outputFolder filesep 'CHN' num2str(channels(h), '%.2d') '.xml'];
        attemptCount = 0;
        while attemptCount ~= -1 && attemptCount < 100
            attemptCount = attemptCount + 1;
            system(sprintf('copy "%s" "%s"', source, target));
            if exist(target, 'file')
                attemptCount = -1;
            else
                pause(0.3);
            end;
        end;
        if attemptCount ~= -1
            error('Failed to copy XML file.');
        end;
    end;
    
    if isempty(startsTop) || isempty(heights) || isempty(startsLeft) || isempty(widths) || isempty(startsFront) || isempty(depths)
        if ~isempty(dimensions)
            stackDimensions = dimensions;
        else
            xmlName = [inputFolder filesep 'ch' num2str(channels(1)) '.xml'];
            attemptCount = 0;
            while attemptCount ~= -1 && attemptCount < 100
                try
                    attemptCount = attemptCount + 1;
                    stackDimensions = readDimensionsFromXML(xmlName);
                    attemptCount = -1;
                catch errorMessage
                    pause(0.3);
                end;
            end;
            if attemptCount ~= -1
                error('Failed to open XML file.');
            end;
            
            if size(stackDimensions, 1) > numel(cameras)
                if size(stackDimensions, 1) >= (max(cameras) + 1)
                    stackDimensions = stackDimensions(cameras + 1, :);
                else
                    error('Unable to retrieve stack dimensions for all cameras from XML files.');
                end;
            elseif size(stackDimensions, 1) < numel(cameras)
                error('Unable to retrieve stack dimensions for all cameras from XML files.');
            end;
        end;
        
        if isempty(startsTop) && isempty(heights)
            startsTop = zeros(1, numel(cameras));
            heights = stackDimensions(:, 2)';
        elseif isempty(startsTop)
            startsTop = zeros(1, numel(cameras));
        elseif isempty(heights)
            heights = stackDimensions(:, 2)' - startsTop;
        end;
        
        if isempty(startsLeft) && isempty(widths)
            startsLeft = zeros(1, numel(cameras));
            widths = stackDimensions(:, 1)';
        elseif isempty(startsLeft)
            startsLeft = zeros(1, numel(cameras));
        elseif isempty(widths)
            widths = stackDimensions(:, 1)' - startsLeft;
        end;
        
        if isempty(startsFront) && isempty(depths)
            startsFront = zeros(1, numel(cameras));
            depths = stackDimensions(:, 3)';
        elseif isempty(startsFront)
            startsFront = zeros(1, numel(cameras));
        elseif isempty(depths)
            depths = stackDimensions(:, 3)' - startsFront;
        end;
    end;
    
    for c = 1:numel(cameras)
        currentCamera = cameras(c);
        
        for h = 1:numel(channels)
            currentChannel = channels(h);
            
            % load stack -> stack
            
            stackName = [inputFolder filesep 'File' num2str(timepoint, '%.4d') ...
                '_CM' num2str(currentCamera) '_CHN' num2str(currentChannel, '%.2d') '.stack'];
            xmlName = [inputFolder filesep 'ch' num2str(currentChannel) '.xml'];
            
            if ~isempty(dimensions)
                stackDimensions = dimensions;
            else
                attemptCount = 0;
                while attemptCount ~= -1 && attemptCount < 100
                    try
                        attemptCount = attemptCount + 1;
                        stackDimensions = readDimensionsFromXML(xmlName);
                        attemptCount = -1;
                    catch errorMessage
                        pause(0.3);
                    end;
                end;
                if attemptCount ~= -1
                    error('Failed to open XML file.');
                end;
                
                if size(stackDimensions, 1) > 1 && (size(stackDimensions, 1) ~= numel(cameras))
                    if size(stackDimensions, 1) >= (currentCamera + 1)
                        stackDimensions = stackDimensions(currentCamera + 1, :);
                    else
                        stackDimensions = stackDimensions(1, :);
                    end;
                elseif size(stackDimensions, 1) > 1
                    stackDimensions = stackDimensions(find(cameras == currentCamera, 1), :);
                end;
            end;
            
            stack = multibandread(stackName, [stackDimensions(2) stackDimensions(1) stackDimensions(3)], ...
                '*uint16', 0, 'bsq', 'ieee-le', ...
                {'Band', 'Range', [startsFront(c) + 1, 1, startsFront(c) + depths(c)]}, ...
                {'Column', 'Range', [startsLeft(c) + 1, 1, startsLeft(c) + widths(c)]}, ...
                {'Row', 'Range', [startsTop(c) + 1, 1, startsTop(c) + heights(c)]});
            
            if rotationFlag == 1
                tempSize = size(stack);
                stack = permute(stack(tempSize(1):-1:1, :, :), [2, 1, 3]);
                clear tempSize;
            elseif rotationFlag == -1
                tempSize = size(stack);
                tempStack = permute(stack, [2, 1, 3]);
                stack = tempStack(tempSize(2):-1:1, :, :);
                clear tempSize tempStack;
            end;
            
            if flipHFlag && c == numel(cameras)
                for z = 1:size(stack, 3)
                    stack(:, :, z) = fliplr(stack(:, :, z));
                end;
            end;
            
            if flipVFlag && c == numel(cameras)
                for z = 1:size(stack, 3)
                    stack(:, :, z) = flipud(stack(:, :, z));
                end;
            end;
            
            % --- MEX file code (does not work yet) ---
            % [fullStack, dimensions] = readBINstack(stackName, int32([height width]));
            % fullStack = reshape(fullStack, dimensions);
            
            % save stack
            
            stackName = [outputFolder filesep 'Sequence' num2str(timepoint, '%.6d') ...
                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentChannel, '%.2d') outputExtension];
            if correctTIFF == 1 && outputType == 2
                writeImage(stack, stackName, 'transpose', 1);
            else
                writeImage(stack, stackName);
            end;
            
            % clear stack
            
            clear stack;
        end;
    end;
    
else % inputType ~= 4
    
    if segmentFlag == 1 || segmentFlag == 2
        kernelSizeArray = [kernelSize  kernelSize  max(1, kernelSize / scaling)];
        kernelSigmaArray = [kernelSigma kernelSigma max(1, kernelSigma / scaling)];
    end;
    
    % generate output sub-directory and copy all XML files to this folder
    
    outputSubFolder = [outputFolder filesep 'TM' num2str(timepoint, '%.6d')];
    if exist(outputSubFolder, 'dir') ~= 7
        mkdir(outputSubFolder);
    end;
    
    if segmentFlag ~= 2
        save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d') '.configuration.mat'], 'configuration');
    end;
    
    if segmentFlag ~= 3
        for h = 1:numel(channels)
            if inputType == 3
                source = [inputFolder filesep 'ch' num2str(channels(h)) '.xml'];
            else
                source = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep 'ch' num2str(channels(h)) '.xml'];
            end;
            target = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d') '_CHN' num2str(channels(h), '%.2d') '.xml'];
            copyfile(source, target);
        end;
    end;
    
    if segmentFlag == 3
        for c = 1:numel(cameras)
            for h = 1:numel(channels)
                % source = [globalMaskFolder filesep 'Global.segmentationMask' outputExtension];
                % target = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') ...
                %     '_CM' num2str(cameras(c), '%.2d') '_CHN' num2str(channels(h), '%.2d') '.segmentationMask' outputExtension];
                % copyfile(source, target);
                
                delete([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                    '_CM' num2str(cameras(c), '%.2d') '_CHN' num2str(channels(h), '%.2d') '.segmentationMask' outputExtension]);
                
                source = [globalMaskFolder filesep 'Global.xzMask' outputExtension];
                target = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                    '_CM' num2str(cameras(c), '%.2d') '_CHN' num2str(channels(h), '%.2d') '.xzMask' outputExtension];
                copyfile(source, target);
                
                source = [globalMaskFolder filesep 'Global.xyMask' outputExtension];
                target = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                    '_CM' num2str(cameras(c), '%.2d') '_CHN' num2str(channels(h), '%.2d') '.xyMask' outputExtension];
                copyfile(source, target);
            end;
        end;
        
        globalMask = readImage([globalMaskFolder filesep 'Global.segmentationMask' outputExtension]);
    end;
    
    if isempty(references)
        references = channels';
        dependents = [];
        if numel(thresholds) ~= numel(references)
            thresholds = ones(numel(references), 1) .* thresholds(1);
        end;
    end;
    
    if isempty(startsTop) || isempty(heights) || isempty(startsLeft) || isempty(widths) || isempty(startsFront) || isempty(depths)
        xmlName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
            '_CHN' num2str(references(1), '%.2d') '.xml'];
        
        if ~isempty(dimensions)
            stackDimensions = dimensions;
        else
            stackDimensions = readDimensionsFromXML(xmlName);
            
            if size(stackDimensions, 1) > numel(cameras)
                if size(stackDimensions, 1) >= (max(cameras) + 1)
                    stackDimensions = stackDimensions(cameras + 1, :);
                else
                    error('Unable to retrieve stack dimensions for all cameras from XML files.');
                end;
            elseif size(stackDimensions, 1) < numel(cameras)
                error('Unable to retrieve stack dimensions for all cameras from XML files.');
            end;
        end;
        
        if isempty(startsTop) && isempty(heights)
            startsTop = zeros(1, numel(cameras));
            heights = stackDimensions(:, 2)';
        elseif isempty(startsTop)
            startsTop = zeros(1, numel(cameras));
        elseif isempty(heights)
            heights = stackDimensions(:, 2)' - startsTop;
        end;
        
        if isempty(startsLeft) && isempty(widths)
            startsLeft = zeros(1, numel(cameras));
            widths = stackDimensions(:, 1)';
        elseif isempty(startsLeft)
            startsLeft = zeros(1, numel(cameras));
        elseif isempty(widths)
            widths = stackDimensions(:, 1)' - startsLeft;
        end;
        
        if isempty(startsFront) && isempty(depths)
            startsFront = zeros(1, numel(cameras));
            depths = stackDimensions(:, 3)';
        elseif isempty(startsFront)
            startsFront = zeros(1, numel(cameras));
        elseif isempty(depths)
            depths = stackDimensions(:, 3)' - startsFront;
        end;
    end;
    
    for c = 1:numel(cameras)
        currentCamera = cameras(c);
        
        for r = 1:size(references, 1)
            
            % determine current group of references, dependents and thresholds
            
            currentReferences = references(r, :);
            nCurrentReferences = numel(currentReferences);
            if ~isempty(dependents)
                currentDependents = dependents(r, :);
            else
                currentDependents = [];
            end;
            nCurrentDependents = numel(currentDependents);
            currentThreshold = thresholds(r);
            
            referenceStacks = cell(nCurrentReferences, 1);
            referenceMasks = cell(nCurrentReferences, 1);
            
            for n = 1:nCurrentReferences
                
                % load stack -> referenceStacks
                
                if inputType == 0 || inputType == 1
                    referenceStacks{n} = zeros(heights(c), widths(c), depths(c), 'uint16');
                    zRange = startsFront(c):(startsFront(c) + depths(c) - 1);
                end;
                
                if inputType == 0
                    for z = 1:numel(zRange)
                        imageName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') '_ANG' num2str(angle, '%.3d') ...
                            '_CM' num2str(currentCamera) '_CHN' num2str(currentReferences(n), '%.2d') '_PH0_PLN' num2str(zRange(z), '%.4d') '.tif'];
                        currentImage = readImage(imageName);
                        referenceStacks{n}(:, :, z) = currentImage((startsTop(c) + 1):(startsTop(c) + heights(c)), ...
                            ((startsLeft(c) + 1):(startsLeft(c) + widths(c))));
                    end;
                elseif inputType == 1
                    for z = 1:numel(zRange)
                        imageName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') '_ANG' num2str(angle, '%.3d') ...
                            '_CM' num2str(currentCamera) '_CHN' num2str(currentReferences(n), '%.2d') '_PH0_PLN' num2str(zRange(z), '%.4d') '.jp2'];
                        currentImage = readImage(imageName);
                        referenceStacks{n}(:, :, z) = currentImage((startsTop(c) + 1):(startsTop(c) + heights(c)), ...
                            ((startsLeft(c) + 1):(startsLeft(c) + widths(c))));
                    end;
                elseif inputType == 2 || inputType == 3
                    if inputType == 2
                        stackName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') '_ANG' num2str(angle, '%.3d') ...
                            '_CM' num2str(currentCamera) '_CHN' num2str(currentReferences(n), '%.2d') '_PH0.stack'];
                    else
                        stackName = [inputFolder filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') '_ANG' num2str(angle, '%.3d') ...
                            '_CM' num2str(currentCamera) '_CHN' num2str(currentReferences(n), '%.2d') '_PH0.stack'];
                    end;
                    xmlName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CHN' num2str(currentReferences(n), '%.2d') '.xml'];
                    
                    if ~isempty(dimensions)
                        stackDimensions = dimensions;
                    else
                        stackDimensions = readDimensionsFromXML(xmlName);
                        if size(stackDimensions, 1) > 1 && (size(stackDimensions, 1) ~= numel(cameras))
                            if size(stackDimensions, 1) >= (currentCamera + 1)
                                stackDimensions = stackDimensions(currentCamera + 1, :);
                            else
                                stackDimensions = stackDimensions(1, :);
                            end;
                        elseif size(stackDimensions, 1) > 1
                            stackDimensions = stackDimensions(find(cameras == currentCamera, 1), :);
                        end;
                    end;
                    
                    referenceStacks{n} = multibandread(stackName, [stackDimensions(2) stackDimensions(1) stackDimensions(3)], ...
                        '*uint16', 0, 'bsq', 'ieee-le', ...
                        {'Band', 'Range', [startsFront(c) + 1, 1, startsFront(c) + depths(c)]}, ...
                        {'Column', 'Range', [startsLeft(c) + 1, 1, startsLeft(c) + widths(c)]}, ...
                        {'Row', 'Range', [startsTop(c) + 1, 1, startsTop(c) + heights(c)]});
                    
                    % --- MEX file code (does not work yet) ---
                    % [fullStack, dimensions] = readBINstack(stackName, int32([height width]));
                    % fullStack = reshape(fullStack, dimensions);
                end;
                
                if rotationFlag == 1
                    tempSize = size(referenceStacks{n});
                    referenceStacks{n} = permute(referenceStacks{n}(tempSize(1):-1:1, :, :), [2, 1, 3]);
                    clear tempSize;
                elseif rotationFlag == -1
                    tempSize = size(referenceStacks{n});
                    tempStack = permute(referenceStacks{n}, [2, 1, 3]);
                    referenceStacks{n} = tempStack(tempSize(2):-1:1, :, :);
                    clear tempSize tempStack;
                end;
                
                if flipHFlag && c == numel(cameras)
                    for z = 1:size(referenceStacks{n}, 3)
                        referenceStacks{n}(:, :, z) = fliplr(referenceStacks{n}(:, :, z));
                    end;
                end;
                
                if flipVFlag && c == numel(cameras)
                    for z = 1:size(referenceStacks{n}, 3)
                        referenceStacks{n}(:, :, z) = flipud(referenceStacks{n}(:, :, z));
                    end;
                end;
                
                % correct stack -> referenceStacks (overwrite)
                
                if ~isempty(medianRange)
                    [referenceStacks{n}, pixelCorrection] = correctInsensitivePixels(referenceStacks{n}, backgroundValues(c), medianRange, verbose);
                end;
                
                if loggingFlag
                    if ~isempty(medianRange)
                        correctionDatabase = {backgroundValues(c), pixelCorrection{1}, pixelCorrection{2}, pixelCorrection{3}, pixelCorrection{4}};
                    else
                        correctionDatabase = backgroundValues(c);
                    end;
                    save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '.correctionDatabase.mat'], 'correctionDatabase');
                    clear correctionDatabase;
                end;
                
                % calculate background [prctile(subSampledStack(subSampledStack > 0), percentile);] and save
                
                if segmentFlag ~= 3
                    subSampledStack = referenceStacks{n}(1:percentile(3):end);
                    if segmentFlag == 1 || segmentFlag == 2
                        minIntensity = zeros(1, 2);
                        minIntensity(2) = prctile(subSampledStack(subSampledStack > 0), percentile(2));
                    else
                        minIntensity = prctile(subSampledStack(subSampledStack > 0), percentile(2));
                        save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                            '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '.minIntensity.mat'], 'minIntensity');
                    end;
                    clear subSampledStack;
                end;
                
                % generate gauss stack -> referenceMasks
                
                if segmentFlag == 1 || segmentFlag == 2
                    if splitting > 1
                        if rotationFlag ~= 0
                            referenceMasks{n} = zeros(widths(c), heights(c), depths(c), 'uint16');
                        else
                            referenceMasks{n} = zeros(heights(c), widths(c), depths(c), 'uint16');
                        end;
                        
                        splittingMargin = 2 * kernelSize;
                        
                        for i = 1:splitting
                            if rotationFlag ~= 0
                                slabStart = max(1, round((i - 1) * widths(c) / splitting + 1 - splittingMargin));
                                slabStop = min(widths(c), round(i * widths(c) / splitting + splittingMargin));
                            else
                                slabStart = max(1, round((i - 1) * heights(c) / splitting + 1 - splittingMargin));
                                slabStop = min(heights(c), round(i * heights(c) / splitting + splittingMargin));
                            end;
                            convolvedSlab = uint16(imgaussianAnisotropy(double(referenceStacks{n}(slabStart:slabStop, :, :)), kernelSigmaArray, kernelSizeArray));
                            if i == 1
                                referenceMasks{n}(1:(slabStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                            elseif i == splitting
                                referenceMasks{n}((slabStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                            else % i > 1 && i < splitting
                                referenceMasks{n}((slabStart + splittingMargin):(slabStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                            end;
                            clear convolvedSlab;
                        end;
                    else
                        referenceMasks{n} = uint16(imgaussianAnisotropy(double(referenceStacks{n}), kernelSigmaArray, kernelSizeArray));
                    end;
                    
                    % generate thresholded mask -> referenceMasks (overwrite)
                    
                    minIntensity(1) = prctile(referenceMasks{n}(1:percentile(3):end), percentile(1));
                    save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d') ...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '.minIntensity.mat'], 'minIntensity');
                    
                    intensityStatistics = zeros(splitting, 2);
                    for i = 1:splitting
                        if rotationFlag ~= 0
                            slabStart = round((i - 1) * widths(c) / splitting + 1);
                            slabStop  = round(i * widths(c) / splitting);
                        else
                            slabStart = round((i - 1) * heights(c) / splitting + 1);
                            slabStop  = round(i * heights(c) / splitting);
                        end;
                        temporaryArray = referenceMasks{n}(slabStart:slabStop, :, :);
                        temporaryArray = temporaryArray(temporaryArray > minIntensity(1));
                        intensityStatistics(i, 1) = sum(temporaryArray(:));
                        intensityStatistics(i, 2) = size(temporaryArray, 1);
                    end;
                    meanIntensity = sum(intensityStatistics(:, 1)) / sum(intensityStatistics(:, 2));
                    clear temporaryArray intensityStatistics;
                    
                    level = minIntensity(1) + (meanIntensity - minIntensity(1)) * currentThreshold;
                    
                    referenceMasks{n} = referenceMasks{n} > level;
                end;
            end;
            
            % if nCurrentReferences > 1 -> fuse referenceMasks using OR -> referenceMasks (overwrite)
            
            currentMasterList = cat(2, currentReferences, currentDependents);
            if nCurrentReferences > 1 && (segmentFlag == 1 || segmentFlag == 2)
                for n = 2:nCurrentReferences
                    referenceMasks{1} = referenceMasks{1} | referenceMasks{n};
                    referenceMasks{n} = [];
                end;
            end;
            
            if segmentFlag == 1 || segmentFlag == 2
                referenceMasks{1} = uint16(referenceMasks{1});
                
                if loggingFlag || segmentFlag == 2
                    segmentationMaskName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(1), '%.2d') '.segmentationMask' outputExtension];
                    if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                        writeImage(referenceMasks{1}, segmentationMaskName, 'transpose', 1);
                    else
                        writeImage(referenceMasks{1}, segmentationMaskName);
                    end;
                    
                    if numel(currentMasterList) > 1 && segmentFlag == 1
                        for n = 2:numel(currentMasterList)
                            targetName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(n), '%.2d') '.segmentationMask' outputExtension];
                            copyfile(segmentationMaskName, targetName);
                        end;
                    end;
                end;
                
                if segmentFlag == 1
                    for n = 1:nCurrentReferences
                        
                        % mask referenceStacks using referenceMasks -> referenceStacks (overwrite)
                        
                        referenceStacks{n} = referenceStacks{n} .* referenceMasks{1};
                        
                        %% Inactive code section for padding option
                        % if padding(1) == 0
                        %     referenceStacks{n} = referenceStacks{n} .* referenceMasks{1};
                        % end;
                    end;
                    
                    % generate coordinate masks using referenceMasks and save (for all references and dependents)
                    
                    if rotationFlag ~= 0
                        xzSliceMask = zeros(widths(c), depths(c), 'uint16');
                    else
                        xzSliceMask = zeros(heights(c), depths(c), 'uint16');
                    end;
                    for i = 1:splitting
                        if rotationFlag ~= 0
                            slabStart = round((i - 1) * widths(c) / splitting + 1);
                            slabStop  = round(i * widths(c) / splitting);
                        else
                            slabStart = round((i - 1) * heights(c) / splitting + 1);
                            slabStop  = round(i * heights(c) / splitting);
                        end;
                        if rotationFlag ~= 0
                            coordinateMask = repmat(1:heights(c), [(slabStop - slabStart + 1) 1 depths(c)]);
                        else
                            coordinateMask = repmat(1:widths(c), [(slabStop - slabStart + 1) 1 depths(c)]);
                        end;
                        coordinateMask(referenceMasks{1}(slabStart:slabStop, :, :) == 0) = NaN;
                        xzSliceMask(slabStart:slabStop, :) = uint16(round(squeeze(nanmean(coordinateMask, 2))));
                    end;
                    clear coordinateMask;
                    
                    %% Inactive code section for padding option
                    % if padding(1) == 1
                    %     if rotationFlag ~= 0
                    %         xzSliceMask(xzSliceMask == 0) = round(heights(c) / 2);
                    %     else
                    %         xzSliceMask(xzSliceMask == 0) = round(widths(c) / 2);
                    %     end;
                    % elseif padding(1) == 2
                    %     % step 1: dilate binary mask by a certain amount (disk element with radius defined by padding(2))
                    %     % step 2: fill region corresponding to outside of dilated mask with center coordinates
                    %     % step 3: interpolate transition zone between original greyscale foreground and new background (center coordinates) using scatteredInterpolant
                    %     if rotationFlag ~= 0
                    %         xzBinaryMask = xzSliceMask > 0;
                    %         xzDilatedMask = imdilate(xzBinaryMask, strel('disk', padding(2)));
                    %         xzSliceMask(~xzDilatedMask) = round(heights(c) / 2);
                    %
                    %         [A, B] = meshgrid(1:depths(c), 1:widths(c));
                    %         foregroundSlots = find(xzSliceMask > 0);
                    %         transitionSlots = find(xzSliceMask == 0);
                    %         f = scatteredInterpolant(A(foregroundSlots), B(foregroundSlots), double(xzSliceMask(foregroundSlots)), 'natural');
                    %         transitionValues = f(A(transitionSlots), B(transitionSlots));
                    %         xzSliceMask(transitionSlots) = transitionValues;
                    %         xzSliceMask(xzSliceMask < 1) = 1;
                    %         xzSliceMask(xzSliceMask > heights(c)) = heights(c);
                    %     else
                    %         xzBinaryMask = xzSliceMask > 0;
                    %         xzDilatedMask = imdilate(xzBinaryMask, strel('disk', padding(2)));
                    %         xzSliceMask(~xzDilatedMask) = round(widths(c) / 2);
                    %
                    %         [A, B] = meshgrid(1:depths(c), 1:heights(c));
                    %         foregroundSlots = find(xzSliceMask > 0);
                    %         transitionSlots = find(xzSliceMask == 0);
                    %         f = scatteredInterpolant(A(foregroundSlots), B(foregroundSlots), double(xzSliceMask(foregroundSlots)), 'natural');
                    %         transitionValues = f(A(transitionSlots), B(transitionSlots));
                    %         xzSliceMask(transitionSlots) = transitionValues;
                    %         xzSliceMask(xzSliceMask < 1) = 1;
                    %         xzSliceMask(xzSliceMask > widths(c)) = widths(c);
                    %     end;
                    % end;
                    
                    if rotationFlag ~= 0
                        xySliceMask = zeros(widths(c), heights(c), 'uint16');
                    else
                        xySliceMask = zeros(heights(c), widths(c), 'uint16');
                    end;
                    for i = 1:splitting
                        if rotationFlag ~= 0
                            slabStart = round((i - 1) * widths(c) / splitting + 1);
                            slabStop  = round(i * widths(c) / splitting);
                        else
                            slabStart = round((i - 1) * heights(c) / splitting + 1);
                            slabStop  = round(i * heights(c) / splitting);
                        end;
                        if rotationFlag ~= 0
                            coordinateMask = repmat(reshape(1:depths(c), [1 1 depths(c)]), [(slabStop - slabStart + 1) heights(c) 1]);
                        else
                            coordinateMask = repmat(reshape(1:depths(c), [1 1 depths(c)]), [(slabStop - slabStart + 1) widths(c) 1]);
                        end;
                        coordinateMask(referenceMasks{1}(slabStart:slabStop, :, :) == 0) = NaN;
                        xySliceMask(slabStart:slabStop, :) = uint16(round(squeeze(nanmean(coordinateMask, 3))));
                    end;
                    clear coordinateMask;
                    
                    %% Inactive code section for padding option
                    % if padding(1) == 1
                    %     xySliceMask(xySliceMask == 0) = round(depths(c) / 2);
                    % elseif padding(1) == 2
                    %     if rotationFlag ~= 0
                    %         xyBinaryMask = xySliceMask > 0;
                    %         xyDilatedMask = imdilate(xyBinaryMask, strel('disk', padding(2)));
                    %         xySliceMask(~xyDilatedMask) = round(depths(c) / 2);
                    %
                    %         [A, B] = meshgrid(1:heights(c), 1:widths(c));
                    %         foregroundSlots = find(xySliceMask > 0);
                    %         transitionSlots = find(xySliceMask == 0);
                    %         f = scatteredInterpolant(A(foregroundSlots), B(foregroundSlots), double(xySliceMask(foregroundSlots)), 'natural');
                    %         transitionValues = f(A(transitionSlots), B(transitionSlots));
                    %         xySliceMask(transitionSlots) = transitionValues;
                    %         xySliceMask(xySliceMask < 1) = 1;
                    %         xySliceMask(xySliceMask > depths(c)) = depths(c);
                    %     else
                    %         xyBinaryMask = xySliceMask > 0;
                    %         xyDilatedMask = imdilate(xyBinaryMask, strel('disk', padding(2)));
                    %         xySliceMask(~xyDilatedMask) = round(depths(c) / 2);
                    %
                    %         [A, B] = meshgrid(1:widths(c), 1:heights(c));
                    %         foregroundSlots = find(xySliceMask > 0);
                    %         transitionSlots = find(xySliceMask == 0);
                    %         f = scatteredInterpolant(A(foregroundSlots), B(foregroundSlots), double(xySliceMask(foregroundSlots)), 'natural');
                    %         transitionValues = f(A(transitionSlots), B(transitionSlots));
                    %         xySliceMask(transitionSlots) = transitionValues;
                    %         xySliceMask(xySliceMask < 1) = 1;
                    %         xySliceMask(xySliceMask > depths(c)) = depths(c);
                    %     end;
                    % end;
                    
                    xzMaskName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(1), '%.2d') '.xzMask' outputExtension];
                    xyMaskName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(1), '%.2d') '.xyMask' outputExtension];
                    if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                        writeImage(xzSliceMask, xzMaskName, 'transpose', 1);
                        writeImage(xySliceMask, xyMaskName, 'transpose', 1);
                    else
                        writeImage(xzSliceMask, xzMaskName);
                        writeImage(xySliceMask, xyMaskName);
                    end;
                    
                    if numel(currentMasterList) > 1
                        for n = 2:numel(currentMasterList)
                            targetXZName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(n), '%.2d') '.xzMask' outputExtension];
                            targetYZName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentMasterList(n), '%.2d') '.xyMask' outputExtension];
                            copyfile(xzMaskName, targetXZName);
                            copyfile(xyMaskName, targetYZName);
                        end;
                    end;
                end;
            elseif segmentFlag == 3
                for n = 1:nCurrentReferences
                    
                    % mask referenceStacks using globalMask -> referenceStacks (overwrite)
                    
                    referenceStacks{n} = referenceStacks{n} .* globalMask;
                end;
            end;
            
            if segmentFlag ~= 2
                for n = 1:nCurrentReferences
                    
                    % save referenceStacks
                    
                    stackName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                        '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') outputExtension];
                    if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                        writeImage(referenceStacks{n}, stackName, 'transpose', 1);
                    else
                        writeImage(referenceStacks{n}, stackName);
                    end;
                    
                    % save 2D projections of referenceStacks
                    for currentFolder = {projectionFolder, [outputFolder filesep 'TM' num2str(timepoint, '%.6d')]}
                        xyProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                            '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '_xyProjection' outputExtension];
                        xzProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                            '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '_xzProjection' outputExtension];
                        yzProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                            '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentReferences(n), '%.2d') '_yzProjection' outputExtension];
                        if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                            writeImage(squeeze(max(referenceStacks{n}, [], 3)), xyProjectionName, 'transpose', 1);
                            writeImage(squeeze(max(referenceStacks{n}, [], 2)), xzProjectionName, 'transpose', 1);
                            writeImage(squeeze(max(referenceStacks{n}, [], 1)), yzProjectionName, 'transpose', 1);
                        else
                            writeImage(squeeze(max(referenceStacks{n}, [], 3)), xyProjectionName);
                            writeImage(squeeze(max(referenceStacks{n}, [], 2)), xzProjectionName);
                            writeImage(squeeze(max(referenceStacks{n}, [], 1)), yzProjectionName);
                        end;
                    end
                end;
                
                % clear referenceStacks
                
                clear referenceStacks;
                
                % if nCurrentDependents > 0 -> sequentially load each stack, correct stack, save background, mask stack, save stack, clear stack
                
                if nCurrentDependents > 0
                    for n = 1:nCurrentDependents
                        
                        % load stack -> dependentStack
                        
                        if inputType == 0 || inputType == 1
                            dependentStack = zeros(heights(c), widths(c), depths(c), 'uint16');
                            zRange = startsFront(c):(startsFront(c) + depths(c) - 1);
                        end;
                        
                        if inputType == 0
                            for z = 1:numel(zRange)
                                imageName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') ...
                                    '_ANG' num2str(angle, '%.3d') '_CM' num2str(currentCamera) '_CHN' num2str(currentDependents(n), '%.2d') '_PH0_PLN' num2str(zRange(z), '%.4d') '.tif'];
                                currentImage = readImage(imageName);
                                dependentStack(:, :, z) = currentImage((startsTop(c) + 1):(startsTop(c) + heights(c)), ...
                                    ((startsLeft(c) + 1):(startsLeft(c) + widths(c))));
                            end;
                        elseif inputType == 1
                            for z = 1:numel(zRange)
                                imageName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') ...
                                    '_ANG' num2str(angle, '%.3d') '_CM' num2str(currentCamera) '_CHN' num2str(currentDependents(n), '%.2d') '_PH0_PLN' num2str(zRange(z), '%.4d') '.jp2'];
                                currentImage = readImage(imageName);
                                dependentStack(:, :, z) = currentImage((startsTop(c) + 1):(startsTop(c) + heights(c)), ...
                                    ((startsLeft(c) + 1):(startsLeft(c) + widths(c))));
                            end;
                        elseif inputType == 2 || inputType == 3
                            if inputType == 2
                                stackName = [inputFolder filesep 'TM' num2str(timepoint, '%.5d') filesep num2str(angle, 'ANG%.3d') filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d') ...
                                    '_ANG' num2str(angle, '%.3d') '_CM' num2str(currentCamera) '_CHN' num2str(currentDependents(n), '%.2d') '_PH0.stack'];
                            else
                                stackName = [inputFolder filesep 'SPC' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.5d')...
                                    '_ANG' num2str(angle, '%.3d') '_CM' num2str(currentCamera) '_CHN' num2str(currentDependents(n), '%.2d') '_PH0.stack'];
                            end;
                            xmlName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') ...
                                    '_ANG' num2str(angle, '%.3d') '_CHN' num2str(currentDependents(n), '%.2d') '.xml'];
                            
                            if ~isempty(dimensions)
                                stackDimensions = dimensions;
                            else
                                stackDimensions = readDimensionsFromXML(xmlName);
                                if size(stackDimensions, 1) > 1 && (size(stackDimensions, 1) ~= numel(cameras))
                                    if size(stackDimensions, 1) >= (currentCamera + 1)
                                        stackDimensions = stackDimensions(currentCamera + 1, :);
                                    else
                                        stackDimensions = stackDimensions(1, :);
                                    end;
                                elseif size(stackDimensions, 1) > 1
                                    stackDimensions = stackDimensions(find(cameras == currentCamera, 1), :);
                                end;
                            end;
                            
                            dependentStack = multibandread(stackName, [stackDimensions(2) stackDimensions(1) stackDimensions(3)], ...
                                '*uint16', 0, 'bsq', 'ieee-le', ...
                                {'Band', 'Range', [startsFront(c) + 1, 1, startsFront(c) + depths(c)]}, ...
                                {'Column', 'Range', [startsLeft(c) + 1, 1, startsLeft(c) + widths(c)]}, ...
                                {'Row', 'Range', [startsTop(c) + 1, 1, startsTop(c) + heights(c)]});
                            
                            % --- MEX file code (does not work yet) ---
                            % [fullStack, dimensions] = readBINstack(stackName, int32([height width]));
                            % fullStack = reshape(fullStack, dimensions);
                        end;
                        
                        if rotationFlag == 1
                            tempSize = size(dependentStack);
                            dependentStack = permute(dependentStack(tempSize(1):-1:1, :, :), [2, 1, 3]);
                            clear tempSize;
                        elseif rotationFlag == -1
                            tempSize = size(dependentStack);
                            tempStack = permute(dependentStack, [2, 1, 3]);
                            dependentStack = tempStack(tempSize(2):-1:1, :, :);
                            clear tempSize tempStack;
                        end;
                        
                        if flipHFlag && c == numel(cameras)
                            for z = 1:size(dependentStack, 3)
                                dependentStack(:, :, z) = fliplr(dependentStack(:, :, z));
                            end;
                        end;
                        
                        if flipVFlag && c == numel(cameras)
                            for z = 1:size(dependentStack, 3)
                                dependentStack(:, :, z) = flipud(dependentStack(:, :, z));
                            end;
                        end;
                        
                        % correct stack -> dependentStack (overwrite)
                        
                        if ~isempty(medianRange)
                            [dependentStack, pixelCorrection] = correctInsensitivePixels(dependentStack, backgroundValues(c), medianRange, verbose);
                        end;
                        
                        if loggingFlag
                            if ~isempty(medianRange)
                                correctionDatabase = {backgroundValues(c), pixelCorrection{1}, pixelCorrection{2}, pixelCorrection{3}, pixelCorrection{4}};
                            else
                                correctionDatabase = backgroundValues(c);
                            end;
                            save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') '.correctionDatabase.mat'], 'correctionDatabase');
                            clear correctionDatabase;
                        end;
                        
                        % calculate background [prctile(subSampledStack(subSampledStack > 0), percentile);] and save
                        
                        if segmentFlag ~= 3
                            subSampledStack = dependentStack(1:percentile(3):end);
                            minIntensity = prctile(subSampledStack(subSampledStack > 0), percentile(2));
                            save([outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') '.minIntensity.mat'], 'minIntensity');
                            clear subSampledStack;
                        end;
                        
                        % mask dependentStack using referenceMasks/globalMask -> dependentStack (overwrite)
                        
                        if segmentFlag == 1
                            dependentStack = dependentStack .* referenceMasks{1};
                        elseif segmentFlag == 3
                            dependentStack = dependentStack .* globalMask;
                        end;
                        
                        %% Inactive code section for padding option
                        % if segmentFlag && padding(1) == 0
                        %     dependentStack = dependentStack .* referenceMasks{1};
                        % end;
                        
                        % save dependentStack
                        
                        stackName = [outputFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                            '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') outputExtension];
                        if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                            writeImage(dependentStack, stackName, 'transpose', 1);
                        else
                            writeImage(dependentStack, stackName);
                        end;
                        
                        % save 2D projections of dependentStack
                        for currentFolder = {projectionFolder, [outputFolder filesep 'TM' num2str(timepoint, '%.6d')]}
                            xyProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d') ...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') '_xyProjection' outputExtension];
                            xzProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') '_xzProjection' outputExtension];
                            yzProjectionName = [currentFolder{:} filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') '_ANG' num2str(angle, '%.3d')...
                                '_CM' num2str(currentCamera, '%.2d') '_CHN' num2str(currentDependents(n), '%.2d') '_yzProjection' outputExtension];
                            if correctTIFF == 1 && ((inputType == 0 && outputType ~= 2) || (inputType ~= 0 && outputType == 2))
                                writeImage(squeeze(max(dependentStack, [], 3)), xyProjectionName, 'transpose', 1);
                                writeImage(squeeze(max(dependentStack, [], 2)), xzProjectionName, 'transpose', 1);
                                writeImage(squeeze(max(dependentStack, [], 1)), yzProjectionName, 'transpose', 1);
                            else
                                writeImage(squeeze(max(dependentStack, [], 3)), xyProjectionName);
                                writeImage(squeeze(max(dependentStack, [], 2)), xzProjectionName);
                                writeImage(squeeze(max(dependentStack, [], 1)), yzProjectionName);
                            end;
                        end;
                        % clear dependentStack
                        
                        clear dependentStack;
                    end;
                end;
            end;
        end;
    end;
end;

end

function dimensions = readDimensionsFromXML(filename)

dimensions = [];

fid = fopen(filename, 'r');
if fid < 1
    error('Could not open file %s for reading.', filename);
end;

while true
    s = fgetl(fid);
    if ~ischar(s)
        break;
    end;
    m = regexp(s, '<info dimensions="([^#]*)"', 'tokens', 'once');
    if ~isempty(m)
        dimensions = m{1};
        break;
    end;
end;
fclose(fid);

if isempty(dimensions)
    error('Unable to retrieve stack dimensions from file %s.', filename);
end

dimensions = str2double(regexp(dimensions, '[^\d]+', 'split'));
dimensions = reshape(dimensions, [3, numel(dimensions) / 3])';

if any(isnan(dimensions))
    error('Unable to correctly parse stack dimensions retrieved from file %s.', filename);
end;

end

function [stack, pixelCorrection] = correctInsensitivePixels(stack, backgroundValue, medianRange, verbose)

deviationProjection = std(single(stack), 0, 3);
deviationProjectionMedianFiltered = medfilt2(deviationProjection, medianRange, 'symmetric');
deviationDistances = abs(deviationProjection - deviationProjectionMedianFiltered);
deviationThreshold = determineThreshold(sort(deviationDistances(:)));

deviationMatrix = deviationDistances > deviationThreshold;

meanProjection = mean(stack, 3) - backgroundValue;
meanProjectionMedianFiltered = medfilt2(meanProjection, medianRange, 'symmetric');
meanDistances = abs((meanProjection - meanProjectionMedianFiltered) ./ meanProjectionMedianFiltered);
meanThreshold = determineThreshold(sort(meanDistances(:)));

meanMatrix = meanDistances > meanThreshold;

pixelMatrix = deviationMatrix | meanMatrix;
pixelCorrection = {deviationDistances, deviationThreshold, meanDistances, meanThreshold};

if verbose
    disp(['Insensitive pixels detected: ' num2str(sum(pixelMatrix(:))) ' (' num2str(100 * sum(pixelMatrix(:)) / numel(pixelMatrix), '%.1f') '%)']);
end;

% apply pixelMatrix to correct insensitive pixels
for z = 1:size(stack, 3)  
    frame = stack(:, :, z);
    filteredFrame = medfilt2(frame, medianRange, 'symmetric');
    frame(pixelMatrix == 1) = filteredFrame(pixelMatrix == 1);
    stack(:, :, z) = frame;
end;

end

function threshold = determineThreshold(array)

warning off; % convergence might not be reached, but solution is still acceptable

elements = length(array);
maxSamples = 50000;

if elements > maxSamples % subsample
    step = round(elements / maxSamples);
    array = array(1:step:end);
    elements = length(array);
end
X = (1:elements)';

% estimate threshold using maximum distance from straight line connecting first and last value
pairs = [X array];
vector = pairs(end, :) - pairs(1, :);
h = dot(repmat(vector, [elements 1]), pairs - repmat(pairs(1, :), [elements 1]), 2) / (norm(vector) ^ 2);
distances = sqrt(sum((repmat(vector, [elements 1]) .* repmat(h, [1 2]) + repmat(pairs(1, :), [elements 1]) - pairs) .^ 2, 2));
[value position] = max(distances);
threshold = array(position);

warning on;

end

function I = imgaussianAnisotropy(I, sigma, siz)

% ---------------------------------------------------------------------------------------------------
% | Anisotropic Gaussian filtering                                                                  |
% | Original function (imgaussian) written by D. Kroon, University of Twente, September 2009        |
% | Code modification to allow different sigmas by Fernando Amat, HHMI/Janelia Farm, September 2010 |
% ---------------------------------------------------------------------------------------------------

% Note: If X is of type double, the MEX file is faster. For images of type single or int, this m-file is faster.

if(~exist('siz', 'var'))
    siz = sigma * 6;
end;

ndimsI = sum(size(I) > 1);
if(length(sigma) ~= ndimsI)
    error 'You must specify one sigma for each dimension of the image'
end;

% Filter each dimension with the 1D Gaussian kernels
if(ndimsI == 1)
    % Make 1D Gaussian kernel
    kk = 1;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));

    I = imfilter(I, H, 'same', 'replicate');
elseif(ndimsI == 2)
    % Make 1D Gaussian kernel
    kk = 1;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));
    Hx = reshape(H, [length(H) 1]);
    
    kk = 2;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));
    Hy = reshape(H, [1 length(H)]);
    
    I = imfilter(imfilter(I, Hx, 'same', 'replicate'), Hy, 'same', 'replicate');
elseif(ndimsI == 3)
    if(size(I, 3) < 4) % Detect if 3D or color image
        % Make 1D Gaussian kernel
        kk = 1;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hx = reshape(H, [length(H) 1]);
        
        kk = 2;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hy = reshape(H, [1 length(H)]);
        for k = 1:size(I, 3)
            I(:, :, k) = imfilter(imfilter(I(:, :, k), Hx, 'same', 'replicate'), Hy, 'same', 'replicate');
        end;
    else
        % Make 1D Gaussian kernel
        kk = 1;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hx = reshape(H, [length(H) 1 1]);
        
        kk = 2;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hy = reshape(H, [1 length(H) 1]);
        
        kk = 3;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hz = reshape(H, [1 1 length(H)]);
        
        I = imfilter(imfilter(imfilter(I, Hx, 'same', 'replicate'), Hy, 'same', 'replicate'), Hz, 'same', 'replicate');
    end;
else
    error('imgaussian:input', 'unsupported input dimension');
end;

end