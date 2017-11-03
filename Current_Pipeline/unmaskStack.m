function unmaskStack(...
    inputFolder, maskedStackName,...
    correctionType, correctionSlot,...
    fastInitialize, erosionRadius, interpRadius, scalingFactor, downsampling,...
    saveMasks, saveMetadata,...
    outputFolder, verbose, nParallelJobs)

% ----------------------------------------------------------------------------
% | Smooth unmasking of image stacks masked by clusterPT/clusterMF/clusterTF |
% |                                                                          |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2016            |
% | Email: kellerp@janelia.hhmi.org                                          |
% ----------------------------------------------------------------------------

% correctionType = 1;              % 1: adaptive unmasking (smooth transition to background intensity)
%                                  % 2: unmasking using fixed background intensity
% correctionSlot = [2 100];        % first value set to 1 or 2: first value identifies minIntensity slot used as background intensity, second value ignored
%                                  % first value set to 3: second value defines background intenstiy (override minIntensity data)

% fastInitialize = 1;              % (correctionType == 1 only) flag for using only surface foreground voxels when executing scatteredInterpolant
% erosionRadius  = 2;              % (correctionType == 1 only) radius of erosion kernel used to find foreground voxels relevant for scatteredInterpolant
% interpRadius   = 5;              % (correctionType == 1 only) size of smooth transition zone bridging foreground and background
% scalingFactor  = 2.031/(6.5/16); % (correctionType == 1 only) ratio of axial-vs-lateral voxel size
% downsampling   = 3;              % (correctionType == 1 only) if >1, this is the downsampling factor applied to image stack for constructing interpolation function
% saveMasks      = 1;              % (correctionType == 1 only) flag for saving intermediate binary and dilated masks
% saveMetadata   = 1;              % (correctionType == 1 only) flag for saving voxel count and timing information

% outputFolder   = 'Unmasked';     % target folder for saving results (unmasked stacks, also intermediate mask results if saveMasks == 1)
% verbose        = 1;              % flag for displaying information on current processing progress
% nParallelJobs  = 4;              % number of jobs that are processed in parallel (set to 1 to disable parallel processing)

% ---
% Notes on software functionality that could be added in the future:
% ---
% 1) Dither image in transition zone to suppress deconvolution artifacts
% 2) Apply median filter to surface voxels before scattered interpolation
% 3) Check for boundary artifacts after interpolation and remove these by post-processing
% ---

if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end;

if fastInitialize
    unmaskedStackName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
        'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.klb'];
else
    unmaskedStackName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
        'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '.klb'];
end;

if exist(unmaskedStackName, 'file') == 2
    disp(['Previous results found for stack ' maskedStackName]);
else
    disp(['Processing stack ' maskedStackName]);
    
    stack = readImage([inputFolder filesep '' maskedStackName]);
    
    if correctionType == 1
        binaryMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) 'binaryMask.klb'];
        dilatedMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) 'dilatedMask_r' num2str(interpRadius) '.klb'];
        
        if exist(dilatedMaskName, 'file') == 2
            if verbose
                disp('*** Loading dilated mask');
            end;
            tic;
            dilatedMask = logical(readImage(dilatedMaskName));
            dilatedMaskTime = toc;
            if verbose
                disp(['    Dilated mask loaded in ' num2str(round(dilatedMaskTime)) ' s']);
            end;
        else
            if verbose
                disp('*** Creating dilated mask');
            end;
            tic;
            binaryMask = stack > 0;
            [xGrid, yGrid, zGrid] = meshgrid(-interpRadius:interpRadius);
            ellipsoid = (sqrt(xGrid.^2 + yGrid.^2 + (zGrid .* scalingFactor).^2) <= interpRadius);
            while ~isempty(ellipsoid) && sum(sum(ellipsoid(:, :, 1))) == 0 && sum(sum(ellipsoid(:, :, end))) == 0
                ellipsoid(:, :, [1 end]) = [];
            end;
            dilatedMask = imdilate(binaryMask, ellipsoid);
            
            if saveMasks
                if exist(binaryMaskName, 'file') ~= 2
                    writeImage(uint8(binaryMask) .* 255, binaryMaskName);
                end;
                writeImage(uint8(dilatedMask) .* 255, dilatedMaskName);
            end;
            if nParallelJobs == 1
                clear binaryMask ellipsoid xGrid yGrid zGrid;
            end;
            dilatedMaskTime = toc;
            if verbose
                disp(['    Dilated mask created in ' num2str(round(dilatedMaskTime)) ' s']);
            end;
        end;
        
        if verbose
            disp('*** Flooding background');
        end;
        tic;
        if correctionSlot(1) == 3
            backgroundLevel = correctionSlot(2);
        else
            load([inputFolder filesep '' maskedStackName(1:(end-3)) 'minIntensity.mat']);
            backgroundLevel = minIntensity(correctionSlot(1));
        end;
        stack(~dilatedMask) = backgroundLevel;
        stack(stack > 0 & stack < backgroundLevel) = backgroundLevel;
        if nParallelJobs == 1
            clear dilatedMask;
        end;
        backgroundFloodingTime = toc;
        if verbose
            disp(['    Background flooded in ' num2str(round(backgroundFloodingTime)) ' s']);
        end;
        
        if downsampling <= 1
            if fastInitialize
                erodedMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
                    'erodedMask_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.klb'];
                minimalForegroundMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
                    'minimalForegroundMask_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.klb'];
                if exist(minimalForegroundMaskName, 'file') == 2
                    if verbose
                        disp('*** Loading minimal foreground mask');
                    end;
                    tic;
                    minimalForegroundMask = logical(readImage(minimalForegroundMaskName));
                    minimalForegroundMaskTime = toc;
                    if verbose
                        disp(['    Minimal foreground mask loaded in ' num2str(round(minimalForegroundMaskTime)) ' s']);
                    end;
                else
                    if verbose
                        disp('*** Creating minimal foreground mask');
                    end;
                    tic;
                    binaryMask = stack > 0;
                    [xGrid, yGrid, zGrid] = meshgrid(-erosionRadius:erosionRadius);
                    ball = (sqrt(xGrid.^2 + yGrid.^2 + zGrid.^2) <= erosionRadius);
                    erodedMask = imerode(binaryMask, ball);
                    minimalForegroundMask = xor(binaryMask, erodedMask);
                    
                    if saveMasks
                        writeImage(uint8(erodedMask) .* 255, erodedMaskName);
                        writeImage(uint8(minimalForegroundMask) .* 255, minimalForegroundMaskName);
                    end;
                    if nParallelJobs == 1
                        clear binaryMask erodedMask ball xGrid yGrid zGrid;
                    end;
                    minimalForegroundMaskTime = toc;
                    if verbose
                        disp(['    Minimal foreground mask created in ' num2str(round(minimalForegroundMaskTime)) ' s']);
                    end;
                end;
                foregroundSlots = find(minimalForegroundMask);
            else
                foregroundSlots = find(stack > 0);
            end;
            nForegroundSlots = numel(foregroundSlots);
            transitionSlots = find(stack == 0);
            nTransitionSlots = numel(transitionSlots);
            
            if verbose
                disp(['*** Creating interpolation function using ' num2str(nForegroundSlots) ' foreground voxels']);
            end;
            tic;
            [A, B, C] = ind2sub(size(stack), foregroundSlots);
            f = scatteredInterpolant(A, B, C, double(stack(foregroundSlots)), 'natural');
            functionTime = toc;
            if verbose
                disp(['    Interpolation function created in ' num2str(round(functionTime)) ' s']);
            end;
            
            if verbose
                disp(['*** Flooding transition zone comprising ' num2str(nTransitionSlots) ' voxels']);
            end;
            tic;
            [A, B, C] = ind2sub(size(stack), transitionSlots);
            transitionValues = f(A, B, C);
            stack(transitionSlots) = transitionValues;
            if nParallelJobs == 1
                clear foregroundSlots transitionSlots transitionValues f A B C;
            end;
            transitionFloodingTime = toc;
            if verbose
                disp(['    Transition zone flooded in ' num2str(round(transitionFloodingTime)) ' s']);
            end;
        else
            if verbose
                disp('*** Down-sampling stack');
            end;
            tic;
            dStack = stack((1+floor(downsampling/2)):downsampling:end, (1+floor(downsampling/2)):downsampling:end, :);
            transitionSlots = find(stack == 0);
            nTransitionSlots = numel(transitionSlots);
            samplingTime = toc;
            if verbose
                disp(['    Stack down-sampled in ' num2str(round(samplingTime)) ' s']);
            end;
            
            if fastInitialize
                erodedMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
                    'erodedMask_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.klb'];
                minimalForegroundMaskName = [outputFolder filesep '' maskedStackName(1:(end-3)) ...
                    'minimalForegroundMask_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.klb'];
                if exist(minimalForegroundMaskName, 'file') == 2
                    if verbose
                        disp('*** Loading minimal foreground mask');
                    end;
                    tic;
                    minimalForegroundMask = logical(readImage(minimalForegroundMaskName));
                    minimalForegroundMaskTime = toc;
                    if verbose
                        disp(['    Minimal foreground mask loaded in ' num2str(round(minimalForegroundMaskTime)) ' s']);
                    end;
                else
                    if verbose
                        disp('*** Creating minimal foreground mask');
                    end;
                    tic;
                    dBinaryMask = dStack > 0;
                    [xGrid, yGrid, zGrid] = meshgrid(-erosionRadius:erosionRadius);
                    ball = (sqrt(xGrid.^2 + yGrid.^2 + zGrid.^2) <= erosionRadius);
                    erodedMask = imerode(dBinaryMask, ball);
                    minimalForegroundMask = xor(dBinaryMask, erodedMask);
                    
                    if saveMasks
                        writeImage(uint8(erodedMask) .* 255, erodedMaskName);
                        writeImage(uint8(minimalForegroundMask) .* 255, minimalForegroundMaskName);
                    end;
                    if nParallelJobs == 1
                        clear dBinaryMask erodedMask ball xGrid yGrid zGrid;
                    end;
                    minimalForegroundMaskTime = toc;
                    if verbose
                        disp(['    Minimal foreground mask created in ' num2str(round(minimalForegroundMaskTime)) ' s']);
                    end;
                end;
                foregroundSlots = find(minimalForegroundMask);
            else
                foregroundSlots = find(dStack > 0);
            end;
            nForegroundSlots = numel(foregroundSlots);
            
            if verbose
                disp(['*** Creating interpolation function using ' num2str(nForegroundSlots) ' foreground voxels']);
            end;
            tic;
            [dA, dB, dC] = ind2sub(size(dStack), foregroundSlots);
            dA = (dA-1).*downsampling + floor(downsampling/2) + 1;
            dB = (dB-1).*downsampling + floor(downsampling/2) + 1;
            f = scatteredInterpolant(dA, dB, dC, double(dStack(foregroundSlots)), 'natural');
            if nParallelJobs == 1
                clear foregroundSlots dStack dA dB dC;
            end;
            functionTime = toc;
            if verbose
                disp(['    Interpolation function created in ' num2str(round(functionTime)) ' s']);
            end;
            
            if verbose
                disp(['*** Flooding transition zone comprising ' num2str(nTransitionSlots) ' voxels']);
            end;
            tic;
            [A, B, C] = ind2sub(size(stack), transitionSlots);
            transitionValues = f(A, B, C);
            stack(transitionSlots) = transitionValues;
            if nParallelJobs == 1
                clear transitionSlots transitionValues f A B C;
            end;
            transitionFloodingTime = toc;
            if verbose
                disp(['    Transition zone flooded in ' num2str(round(transitionFloodingTime)) ' s']);
            end;
        end;
    elseif correctionType == 2
        stack(stack == 0) = minIntensity(correctionSlot);
    end;
    writeImage(stack, unmaskedStackName);
    
    if saveMetadata && correctionType == 1
        if downsampling <= 1
            samplingTime = 0;
        end;
        if fastInitialize == 0
            minimalForegroundMaskTime = 0;
        end;
        timeArray = [dilatedMaskTime; backgroundFloodingTime; samplingTime; minimalForegroundMaskTime; functionTime; transitionFloodingTime];
        if fastInitialize
            save([outputFolder filesep '' maskedStackName(1:(end-3)) 'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.timeArray.mat'], 'timeArray');
        else
            save([outputFolder filesep '' maskedStackName(1:(end-3)) 'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '.timeArray.mat'], 'timeArray');
        end;
        
        voxelArray = [nForegroundSlots; nTransitionSlots];
        if fastInitialize
            save([outputFolder filesep '' maskedStackName(1:(end-3)) 'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.voxelArray.mat'], 'voxelArray');
        else
            save([outputFolder filesep '' maskedStackName(1:(end-3)) 'unmasked_r' num2str(interpRadius) '_d' num2str(downsampling) '.voxelArray.mat'], 'voxelArray');
        end;
    end;
end;