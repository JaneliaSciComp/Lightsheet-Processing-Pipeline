%% parameters

inputRoot        = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/Paper/Supplementary_Data_2/Image_Data.MultiFused/';%'X:\SiMView2\13-12-30\Pha_E1_H2BRFP_01_20131230_140802.corrected\Results\MultiFused';
inputPattern     = 'Dme_E1_H2ARFP.TM??????_multiFused_blending/SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';%'Pha_E1_H2BRFP.TM??????_multiFused_blending\SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';
configRoot       = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/Paper/Supplementary_Data_2/Scripts/SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';%'X:\SiMView2\13-12-30\Pha_E1_H2BRFP_01_20131230_140802.corrected\Scripts\SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';

timepoints       = 0:10;%733;
gamma            = 1;
percentile       = 1;

inputType        = 0;       % 0: input data in KLB format
                            % 1: input data in JP2 format
                            % 2: input data in TIF format
outputType       = 0;       % 0: output data saved in KLB format
                            % 1: output data saved in JP2 format
                            % 2: output data saved in TIF format

% configuration of intensity normalization
intensityFlag    = 1;       % flag for enabling/disabling intensity normalization
useStacks        = [1 10];  % 0: use projections for intensity estimate, 1: use stacks for intensity estimate (second parameter provides sub-sampling rate
histogramBins    = 0:(2^16 - 1);
threshold        = 10;      % threshold for histogram computation
backgroundSlot   = 0;       % 0: no background correction, 1: top right, 2: bottom right, 3: bottom left, 4: top left
backgroundEdge   = 100;     % edge size of the square used to estimate background levels from the image data
backgroundDist   = 0;       % vertical/horizontal distance of background square from image reference point defined by backgroundSlot
                            % note: the background square should be positioned in a part of the image without significant foreground content
                            %       throughout the time-lapse experiment

% configuration of correlation-based drift correction
correlationFlag  = 1;       % flag for enabling/disabling correlation-based drift correction

% configuration of global drift correction
globalMode       = 1;       % 0: disable global drift correction
                            % 1: automatic global drift correction based on geometrical center computation
                            % 2: manual global correction using vectors provided for reference time points

% parameters required when globalMode == 1
maskFactor       = 0.4;     % fraction of (mean - minimum) intensity that is used to create the thresholded image mask
scaling          = 2.031 / (6.5 / 16); % axial step size <divided by> (pixel pitch <divided by> magnification)
smoothing        = [1 20];  % smoothing flag ('rloess'), smoothing window size
kernelSize       = 51;      % size of kernel used for Gaussian smoothing prior to thresholding
kernelSigma      = 20;      % sigma of Gaussian smoothing kernel
maskMinimum      = 1;       % percentile for minimum calculation (0 = true minimum)
fraction         = 0.01;    % minimal object size in binary mask, set to 0 to disable bwareaopen (default value 10 ^ -5)

% parameters required when globalMode == 2
referenceDrift   = [...     % each row follows this structure: reference time point (slot 1), x-/y-/z-drift vector (slots 2-4), follows Matlab convention
      0,  0,  0,  0; ...    % Note: at least two rows are required to enable manual global drift correction
    100, 10, 20, 30];

maxStampDigits   = 6;
poolWorkers      = 0;       % use "0" to enable automated detection of available CPU cores
%% from keller
% % % %% parameters
% % % 
% % % inputRoot        = 'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Results' filesep 'MultiFused';
% % % inputPattern     = 'Mmu_E1_CAGTAG1.TM??????_multiFused_blending' filesep 'SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';
% % % configRoot       = 'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Scripts' filesep 'SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';
% % % 
% % % timepoints       = 0:570;
% % % gamma            = 1;
% % % percentile       = 1;
% % % 
% % % inputType        = 0;       % 0: input data in KLB format
% % %                             % 1: input data in JP2 format
% % %                             % 2: input data in TIF format
% % % outputType       = 0;       % 0: output data saved in KLB format
% % %                             % 1: output data saved in JP2 format
% % %                             % 2: output data saved in TIF format
% % % 
% % % % configuration of intensity normalization
% % % intensityFlag    = 1;       % flag for enabling/disabling intensity normalization
% % % useStacks        = [1 10];  % 0: use projections for intensity estimate, 1: use stacks for intensity estimate (second parameter provides sub-sampling rate
% % % histogramBins    = 0:(2^16 - 1);
% % % threshold        = 10;      % threshold for histogram computation
% % % backgroundSlot   = 0;       % 0: no background correction, 1: top right, 2: bottom right, 3: bottom left, 4: top left
% % % backgroundEdge   = 100;
% % % backgroundDist   = 0;
% % % 
% % % % configuration of correlation-based drift correction
% % % correlationFlag  = 1;       % flag for enabling/disabling correlation-based drift correction
% % % 
% % % % configuration of global drift correction
% % % globalMode       = 1;       % 0: disable global drift correction
% % %                             % 1: automatic global drift correction based on geometrical center computation
% % %                             % 2: manual global correction using vectors provided for reference time points
% % % 
% % % % parameters required when globalMode == 1
% % % maskFactor       = 0.4;
% % % scaling          = 2.031 / (6.5 / 16);
% % % smoothing        = [1 20];  % smoothing flag ('rloess'), smoothing window size
% % % kernelSize       = 51;
% % % kernelSigma      = 20;
% % % maskMinimum      = 1;       % percentile for minimum calculation (0 = true minimum)
% % % fraction         = 0.01;    % minimal object size in binary mask, set to 0 to disable bwareaopen (default value 10 ^ -5)
% % % 
% % % % parameters required when globalMode == 2
% % % referenceDrift   = [...     % each row follows this structure: reference time point (slot 1), x-/y-/z-drift vector (slots 2-4), follows Matlab convention
% % %       0,  0,  0,  0; ...    % Note: at least two rows are required to enable manual global drift correction
% % %     100, 10, 20, 30];
% % % 
% % % maxStampDigits   = 6;
% % % poolWorkers      = 0;       % use "0" to enable automated detection of available CPU cores

%% main loop

if poolWorkers == 0
    poolWorkers = feature('numcores');
    disp(' ');
    disp([num2str(poolWorkers) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

if ~exist(configRoot, 'dir')
    mkdir(configRoot);
end;

nTimepoints = numel(timepoints);
inputDatabase = cell(nTimepoints, 1);

for n = 1:nTimepoints
    inputDatabase{n} = [inputRoot filesep '' inputPattern];
    for i = maxStampDigits:-1:1
        positions = strfind(inputDatabase{n}, repmat('?', [1 i]));
        precision = ['%.' num2str(i) 'd'];
        for k = 1:length(positions)
            inputDatabase{n}(positions(k):(positions(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
    end;
end;

disp(' ');
if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);

disp(' ');
disp('Collecting information on image volume size as a function of time');

dimensions = zeros(nTimepoints, 3);

parfor n = 1:nTimepoints
    xyProjectionName = [inputDatabase{n} '_xyProjection' inputExtension];
    xzProjectionName = [inputDatabase{n} '_xzProjection' inputExtension];
    
    switch inputType
        case 0
            headerInformationXY = readKLBheader(xyProjectionName);
            stackDimensionsXY = headerInformationXY.xyzct(1:3);
            
            headerInformationXZ = readKLBheader(xzProjectionName);
            stackDimensionsXZ = headerInformationXZ.xyzct(1:3);
        case 1
            [stackDimensionsXY, bitDepthXY] = readJP2header(xyProjectionName);
            
            [stackDimensionsXZ, bitDepthXZ] = readJP2header(xzProjectionName);
        case 2
            headerInformationXY = imfinfo(xyProjectionName);
            stackDimensionsXY = [headerInformationXY(1).Height headerInformationXY(1).Width numel(headerInformationXY)];
            
            headerInformationXZ = imfinfo(xzProjectionName);
            stackDimensionsXZ = [headerInformationXZ(1).Height headerInformationXZ(1).Width numel(headerInformationXZ)];
    end;
    
    dimensions(n, :) = [stackDimensionsXY(1) stackDimensionsXY(2) stackDimensionsXZ(2)];
end;

dimensions = cat(2, timepoints', dimensions);
save([configRoot filesep 'dimensions.mat'], 'dimensions');

dimensionsMax = max(dimensions(:, 2:4), [], 1);
save([configRoot filesep 'dimensionsMax.mat'], 'dimensionsMax');

if correlationFlag == 1 || globalMode ~= 0
    dimensionsDeltas = zeros(numel(timepoints), 4);
    dimensionsDeltas(:, 1) = timepoints;
    
    for n = 1:nTimepoints
        dimensionsDeltas(n, 2) = floor((dimensionsMax(1) - dimensions(n, 2)) / 2);
        dimensionsDeltas(n, 3) = floor((dimensionsMax(2) - dimensions(n, 3)) / 2);
        dimensionsDeltas(n, 4) = floor((dimensionsMax(3) - dimensions(n, 4)) / 2);
    end;
    
    save([configRoot filesep 'dimensionsDeltas.mat'], 'dimensionsDeltas');
end;

if intensityFlag
    intensityBackgrounds = zeros(nTimepoints, 1);
    intensityCenters = zeros(nTimepoints, 1);
    intensityFactors = zeros(nTimepoints, 1);
    
    save([configRoot filesep 'intensityThreshold.mat'], 'threshold');
end;

if globalMode == 1
    maskStackXY = zeros(dimensionsMax(1), dimensionsMax(2), nTimepoints, 'uint16');
    maskStackXZ = zeros(dimensionsMax(1), dimensionsMax(3), nTimepoints, 'uint16');
    maskStackYZ = zeros(dimensionsMax(2), dimensionsMax(3), nTimepoints, 'uint16');
    
    gaussKernelXY = fspecial('gaussian', [kernelSize, kernelSize], kernelSigma);
    gaussKernelXZ = imresize(gaussKernelXY, [kernelSize, round(kernelSize / scaling)]);
    gaussKernelYZ = gaussKernelXZ;
end;

if correlationFlag == 1 && gamma ~= 1
    gammaStackXY = zeros(dimensionsMax(1), dimensionsMax(2), nTimepoints, 'uint16');
    gammaStackXZ = zeros(dimensionsMax(1), dimensionsMax(3), nTimepoints, 'uint16');
    gammaStackYZ = zeros(dimensionsMax(2), dimensionsMax(3), nTimepoints, 'uint16');
end;

if correlationFlag == 1 || globalMode ~= 0
    driftOffsets = zeros(numel(timepoints), 3);
    driftCenters = zeros(numel(timepoints), 3);
end;

warning off;

if intensityFlag == 1 || correlationFlag == 1 || globalMode == 1
    disp(' ');
    
    parfor n = 1:nTimepoints
        if intensityFlag == 1
            disp(['Estimating intensity for time point ' num2str(timepoints(n), '%.6d')]);
            
            if useStacks(1) == 0
                currentImageXY = readImage([inputDatabase{n} '_xyProjection' inputExtension]);
                
                switch backgroundSlot
                    case 0
                        currentRegion = [];
                    case 1
                        currentRegion = currentImageXY((1 + backgroundDist):(1 + backgroundDist + backgroundEdge), (end - backgroundDist - backgroundEdge + 1):(end - backgroundDist));
                    case 2
                        currentRegion = currentImageXY((end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), (end - backgroundDist - backgroundEdge + 1):(end - backgroundDist));
                    case 3
                        currentRegion = currentImageXY((end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), (1 + backgroundDist):(1 + backgroundDist + backgroundEdge));
                    case 4
                        currentRegion = currentImageXY((1 + backgroundDist):(1 + backgroundDist + backgroundEdge), (1 + backgroundDist):(1 + backgroundDist + backgroundEdge));
                end;
                
                if ~isempty(currentRegion)
                    intensityBackgrounds(n) = mean(currentRegion(:));
                end;
                
                currentHistogram = hist(currentImageXY(:) .* uint16(currentImageXY(:) > (threshold + intensityBackgrounds(n))), histogramBins);
            else
                currentStack = readImage([inputDatabase{n} inputExtension]);
                
                switch backgroundSlot
                    case 0
                        currentRegion = [];
                    case 1
                        currentRegion = currentStack((1 + backgroundDist):(1 + backgroundDist + backgroundEdge), (end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), :);
                    case 2
                        currentRegion = currentStack((end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), (end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), :);
                    case 3
                        currentRegion = currentStack((end - backgroundDist - backgroundEdge + 1):(end - backgroundDist), (1 + backgroundDist):(1 + backgroundDist + backgroundEdge), :);
                    case 4
                        currentRegion = currentStack((1 + backgroundDist):(1 + backgroundDist + backgroundEdge), (1 + backgroundDist):(1 + backgroundDist + backgroundEdge), :);
                end;
                
                if ~isempty(currentRegion)
                    intensityBackgrounds(n) = mean(currentRegion(:));
                end;
                
                currentHistogram = hist(currentStack(1:useStacks(2):end) .* uint16(currentStack(1:useStacks(2):end) > (threshold + intensityBackgrounds(n))), histogramBins);
            end;
            
            currentHistogramCenter = sum(currentHistogram((threshold + 1):end) .* histogramBins((threshold + 1):end)) / sum(currentHistogram((threshold + 1):end));
            intensityCenters(n) = currentHistogramCenter;
        end;
        
        if correlationFlag == 1 || globalMode == 1
            disp(['Estimating drift offsets for time point ' num2str(timepoints(n), '%.6d')]);
            
            if (dimensionsMax(1) ~= dimensions(n, 2)) || (dimensionsMax(2) ~= dimensions(n, 3)) || (dimensionsMax(3) ~= dimensions(n, 4))
                if intensityFlag == 1 && useStacks(1) == 0
                    currentImageXYRaw = currentImageXY;
                else
                    currentImageXYRaw = readImage([inputDatabase{n} '_xyProjection' inputExtension]);
                end;
                currentImageXZRaw = readImage([inputDatabase{n} '_xzProjection' inputExtension]);
                currentImageYZRaw = readImage([inputDatabase{n} '_yzProjection' inputExtension]);
                
                if percentile == 0
                    currentImageXYMin = min(currentImageXYRaw(:));
                    currentImageXZMin = min(currentImageXZRaw(:));
                    currentImageYZMin = min(currentImageYZRaw(:));
                else
                    currentImageXYMin = prctile(currentImageXYRaw(:), percentile);
                    currentImageXZMin = prctile(currentImageXZRaw(:), percentile);
                    currentImageYZMin = prctile(currentImageYZRaw(:), percentile);
                end;
                
                currentImageXY = ones(dimensionsMax(1), dimensionsMax(2), 'uint16') .* currentImageXYMin;
                currentImageXZ = ones(dimensionsMax(1), dimensionsMax(3), 'uint16') .* currentImageXZMin;
                currentImageYZ = ones(dimensionsMax(2), dimensionsMax(3), 'uint16') .* currentImageYZMin;
                
                currentImageXY((dimensionsDeltas(n, 2) + 1):(dimensionsDeltas(n, 2) + dimensions(n, 2)), (dimensionsDeltas(n, 3) + 1):(dimensionsDeltas(n, 3) + dimensions(n, 3))) = currentImageXYRaw;
                currentImageXZ((dimensionsDeltas(n, 2) + 1):(dimensionsDeltas(n, 2) + dimensions(n, 2)), (dimensionsDeltas(n, 4) + 1):(dimensionsDeltas(n, 4) + dimensions(n, 4))) = currentImageXZRaw;
                currentImageYZ((dimensionsDeltas(n, 3) + 1):(dimensionsDeltas(n, 3) + dimensions(n, 3)), (dimensionsDeltas(n, 4) + 1):(dimensionsDeltas(n, 4) + dimensions(n, 4))) = currentImageYZRaw;
            else
                if intensityFlag ~= 1 || useStacks(1) == 1
                    currentImageXY = readImage([inputDatabase{n} '_xyProjection' inputExtension]);
                end;
                currentImageXZ = readImage([inputDatabase{n} '_xzProjection' inputExtension]);
                currentImageYZ = readImage([inputDatabase{n} '_yzProjection' inputExtension]);
            end;
            
            if globalMode == 1 % geometrical center estimation
                filteredImageXY = imfilter(currentImageXY, gaussKernelXY);
                filteredImageXZ = imfilter(currentImageXZ, gaussKernelXZ);
                filteredImageYZ = imfilter(currentImageYZ, gaussKernelYZ);
                
                [xCenter1, yCenter1, xyMask] = calculateCenter(filteredImageXY, maskFactor, maskMinimum, fraction);
                [xCenter2, zCenter1, xzMask] = calculateCenter(filteredImageXZ, maskFactor, maskMinimum, fraction);
                [yCenter2, zCenter2, yzMask] = calculateCenter(filteredImageYZ, maskFactor, maskMinimum, fraction);
                
                maskStackXY(:, :, n) = uint16(xyMask) .* filteredImageXY;
                maskStackXZ(:, :, n) = uint16(xzMask) .* filteredImageXZ;
                maskStackYZ(:, :, n) = uint16(yzMask) .* filteredImageYZ;
                
                xCenter = xCenter1;
                yCenter = yCenter1;
                zCenter = (zCenter1 + zCenter2) / 2;
                driftCenters(n, :) = [xCenter, yCenter, zCenter];
            end;
            
            if correlationFlag == 1 && n > 1 % pairwise correlation-based drift estimation
                if (dimensionsMax(1) ~= dimensions(n - 1, 2)) || (dimensionsMax(2) ~= dimensions(n - 1, 3)) || (dimensionsMax(3) ~= dimensions(n - 1, 4))
                    lastImageXYRaw = readImage([inputDatabase{n - 1} '_xyProjection' inputExtension]);
                    lastImageXZRaw = readImage([inputDatabase{n - 1} '_xzProjection' inputExtension]);
                    lastImageYZRaw = readImage([inputDatabase{n - 1} '_yzProjection' inputExtension]);
                    
                    if percentile == 0
                        lastImageXYMin = min(lastImageXYRaw(:));
                        lastImageXZMin = min(lastImageXZRaw(:));
                        lastImageYZMin = min(lastImageYZRaw(:));
                    else
                        lastImageXYMin = prctile(lastImageXYRaw(:), percentile);
                        lastImageXZMin = prctile(lastImageXZRaw(:), percentile);
                        lastImageYZMin = prctile(lastImageYZRaw(:), percentile);
                    end;
                    
                    lastImageXY = ones(dimensionsMax(1), dimensionsMax(2), 'uint16') .* lastImageXYMin;
                    lastImageXZ = ones(dimensionsMax(1), dimensionsMax(3), 'uint16') .* lastImageXZMin;
                    lastImageYZ = ones(dimensionsMax(2), dimensionsMax(3), 'uint16') .* lastImageYZMin;
                    
                    lastImageXY((dimensionsDeltas(n - 1, 2) + 1):(dimensionsDeltas(n - 1, 2) + dimensions(n - 1, 2)), (dimensionsDeltas(n - 1, 3) + 1):(dimensionsDeltas(n - 1, 3) + dimensions(n - 1, 3))) = lastImageXYRaw;
                    lastImageXZ((dimensionsDeltas(n - 1, 2) + 1):(dimensionsDeltas(n - 1, 2) + dimensions(n - 1, 2)), (dimensionsDeltas(n - 1, 4) + 1):(dimensionsDeltas(n - 1, 4) + dimensions(n - 1, 4))) = lastImageXZRaw;
                    lastImageYZ((dimensionsDeltas(n - 1, 3) + 1):(dimensionsDeltas(n - 1, 3) + dimensions(n - 1, 3)), (dimensionsDeltas(n - 1, 4) + 1):(dimensionsDeltas(n - 1, 4) + dimensions(n - 1, 4))) = lastImageYZRaw;
                else
                    lastImageXY = readImage([inputDatabase{n - 1} '_xyProjection' inputExtension]);
                    lastImageXZ = readImage([inputDatabase{n - 1} '_xzProjection' inputExtension]);
                    lastImageYZ = readImage([inputDatabase{n - 1} '_yzProjection' inputExtension]);
                end;
                
                if gamma ~= 1
                    lastImageXY = imadjust(lastImageXY, [0; 1], [0; 1], gamma);
                    lastImageXZ = imadjust(lastImageXZ, [0; 1], [0; 1], gamma);
                    lastImageYZ = imadjust(lastImageYZ, [0; 1], [0; 1], gamma);
                    
                    currentImageXY = imadjust(currentImageXY, [0; 1], [0; 1], gamma);
                    currentImageXZ = imadjust(currentImageXZ, [0; 1], [0; 1], gamma);
                    currentImageYZ = imadjust(currentImageYZ, [0; 1], [0; 1], gamma);
                    
                    gammaStackXY(:, :, n) = currentImageXY;
                    gammaStackXZ(:, :, n) = currentImageXZ;
                    gammaStackYZ(:, :, n) = currentImageYZ;
                end;
                
                [xOffset1, yOffset1] = correlatePhases(lastImageXY, currentImageXY);
                [xOffset2, zOffset1] = correlatePhases(lastImageXZ, currentImageXZ);
                [yOffset2, zOffset2] = correlatePhases(lastImageYZ, currentImageYZ);
                
                xOffset = xOffset1;
                yOffset = yOffset1;
                zOffset = (zOffset1 + zOffset2) / 2;
                driftOffsets(n, :) = [xOffset, yOffset, zOffset];
            end;
        end;
    end;
    
    disp(' ');
end;

warning on;

if matlabpool('size') > 0
    matlabpool('close');
end;

disp(' ');

if intensityFlag == 1
    maximumSlot = find(intensityCenters == max(intensityCenters), 1);
    referenceHistogramCenter = intensityCenters(maximumSlot);
    intensityFactors = referenceHistogramCenter ./ intensityCenters;
    
    save([configRoot filesep 'intensityBackgrounds.mat'], 'intensityBackgrounds');
    save([configRoot filesep 'intensityCenters.mat'], 'intensityCenters');
    save([configRoot filesep 'intensityFactors.mat'], 'intensityFactors');
    
    figure;
    plot(timepoints, intensityBackgrounds);
    title('Background intensities');
    xlabel('Time point');
    ylabel('Intensity (grey levels)');
    h = gca;
    currentFile = [configRoot filesep 'intensityStep1_intensityBackgrounds.png'];
    saveas(h, currentFile);
    currentFrame = imread(currentFile);
    imwrite(currentFrame, currentFile, 'png');
    
    figure;
    plot(timepoints, intensityCenters);
    title('Intensity histogram centers');
    xlabel('Time point');
    ylabel('Intensity (grey levels)');
    h = gca;
    currentFile = [configRoot filesep 'intensityStep2_intensityCenters.png'];
    saveas(h, currentFile);
    currentFrame = imread(currentFile);
    imwrite(currentFrame, currentFile, 'png');
    
    figure;
    plot(timepoints, intensityFactors);
    title('Intensity correction factors');
    xlabel('Time point');
    ylabel('Factor');
    h = gca;
    currentFile = [configRoot filesep 'intensityStep3_intensityFactors.png'];
    saveas(h, currentFile);
    currentFrame = imread(currentFile);
    imwrite(currentFrame, currentFile, 'png');
end;

if correlationFlag == 1 || globalMode ~= 0
    if correlationFlag == 1 && gamma ~= 1
        disp('Saving gamma-corrected projections');
        disp(' ');
        
        writeImage(gammaStackXY, [configRoot filesep 'xyGamma' outputExtension]);
        writeImage(gammaStackXZ, [configRoot filesep 'xzGamma' outputExtension]);
        writeImage(gammaStackYZ, [configRoot filesep 'yzGamma' outputExtension]);
    end;
    
    cumulativeDriftOffsets = cumsum(driftOffsets, 1);
    
    switch globalMode
        case 0
            smoothedDriftError = zeros(nTimepoints, 3);
        case 1
            disp('Saving masked and blurred projections');
            disp(' ');
            
            writeImage(maskStackXY, [configRoot filesep 'xyMasks' outputExtension]);
            writeImage(maskStackXZ, [configRoot filesep 'xzMasks' outputExtension]);
            writeImage(maskStackYZ, [configRoot filesep 'yzMasks' outputExtension]);
            
            driftCenters = repmat(driftCenters(1, :), nTimepoints, 1) - driftCenters;
            driftFluctuations = cumulativeDriftOffsets - driftCenters;
            
            if smoothing(1) > 0
                smoothedDriftErrorX = smooth(driftFluctuations(:, 1), smoothing(2), 'rloess');
                smoothedDriftErrorY = smooth(driftFluctuations(:, 2), smoothing(2), 'rloess');
                smoothedDriftErrorZ = smooth(driftFluctuations(:, 3), smoothing(2), 'rloess');
                smoothedDriftError = [smoothedDriftErrorX, smoothedDriftErrorY, smoothedDriftErrorZ];
            else
                smoothedDriftError = driftFluctuations;
            end;
        case 2
            nReferenceTimepoints = size(referenceDrift, 1);
            referenceDrift(:, 2:4) = referenceDrift(:, 2:4) - repmat(referenceDrift(1, 2:4), nReferenceTimepoints, 1);
            
            referenceIndices = zeros(nReferenceTimepoints, 1);
            for n = 1:nReferenceTimepoints
                referenceIndices(n) = find(timepoints == referenceDrift(n, 1));
            end;
            
            smoothedDriftErrorX = interp1(referenceDrift(:, 1), cumulativeDriftOffsets(referenceIndices, 1) - referenceDrift(:, 2), timepoints);
            smoothedDriftErrorY = interp1(referenceDrift(:, 1), cumulativeDriftOffsets(referenceIndices, 2) - referenceDrift(:, 3), timepoints);
            smoothedDriftErrorZ = interp1(referenceDrift(:, 1), cumulativeDriftOffsets(referenceIndices, 3) - referenceDrift(:, 4), timepoints);
            smoothedDriftError = [smoothedDriftErrorX; smoothedDriftErrorY; smoothedDriftErrorZ];
            smoothedDriftError = smoothedDriftError';
    end;
    
    finalDriftOffsets = cumulativeDriftOffsets - smoothedDriftError;
    
    driftTable = zeros(nTimepoints, 4);
    driftTable(:, 1) = timepoints;
    driftTable(:, 2:4) = finalDriftOffsets;
    save([configRoot filesep 'driftTable.mat'], 'driftTable');
    
    save([configRoot filesep 'driftDatabase.mat'], 'finalDriftOffsets');
    
    if correlationFlag == 1
        figure;
        plot(timepoints, driftOffsets);
        title('Correlation-based pairwise offsets');
        legend('xOffset', 'yOffset', 'zOffset');
        xlabel('Time point');
        ylabel('Offset (pixel)');
        h = gca;
        currentFile = [configRoot filesep 'driftStep1_pairwiseOffsets.png'];
        saveas(h, currentFile);
        currentFrame = imread(currentFile);
        imwrite(currentFrame, currentFile, 'png');
        
        figure;
        plot(timepoints, cumulativeDriftOffsets);
        title('Cumulative correlation-based offsets');
        legend('xOffset', 'yOffset', 'zOffset');
        xlabel('Time point');
        ylabel('Offset (pixel)');
        h = gca;
        currentFile = [configRoot filesep 'driftStep2_pairwiseOffsetsCumulative.png'];
        saveas(h, currentFile);
        currentFrame = imread(currentFile);
        imwrite(currentFrame, currentFile, 'png');
        
        save([configRoot filesep 'driftDatabase.mat'], 'driftOffsets', 'cumulativeDriftOffsets', '-append');
    end;
    
    if globalMode == 1
        figure;
        plot(timepoints, driftCenters);
        title('Relative geometrical specimen centers');
        legend('xCenter', 'yCenter', 'zCenter');
        xlabel('Time point');
        ylabel('Center (pixel)');
        h = gca;
        currentFile = [configRoot filesep 'driftStep3_specimenCenters.png'];
        saveas(h, currentFile);
        currentFrame = imread(currentFile);
        imwrite(currentFrame, currentFile, 'png');
        
        save([configRoot filesep 'driftDatabase.mat'], 'driftCenters', '-append');
        
        if correlationFlag == 1
            figure;
            plot(timepoints, driftFluctuations);
            title('Drift fluctuations');
            legend('xDrift', 'yDrift', 'zDrift');
            xlabel('Time point');
            ylabel('Drift (pixel)');
            h = gca;
            currentFile = [configRoot filesep 'driftStep4_driftFluctuations.png'];
            saveas(h, currentFile);
            currentFrame = imread(currentFile);
            imwrite(currentFrame, currentFile, 'png');
            
            figure;
            plot(timepoints, smoothedDriftError);
            title('Smoothed drift fluctuations (offset drift)');
            legend('xDrift', 'yDrift', 'zDrift');
            xlabel('Time point');
            ylabel('Drift (pixel)');
            h = gca;
            currentFile = [configRoot filesep 'driftStep5_driftFluctuationsFiltered.png'];
            saveas(h, currentFile);
            currentFrame = imread(currentFile);
            imwrite(currentFrame, currentFile, 'png');
            
            save([configRoot filesep 'driftDatabase.mat'], 'driftFluctuations', 'smoothedDriftError', '-append');
        end;
    elseif globalMode == 2
        figure;
        plot(timepoints, smoothedDriftError);
        if correlationFlag == 0
            title('Relative manual drift correction');
        else
            title('Smoothed drift fluctuations (offset drift)');
        end
        legend('xDrift', 'yDrift', 'zDrift');
        xlabel('Time point');
        ylabel('Drift (pixel)');
        h = gca;
        currentFile = [configRoot filesep 'driftStep3_manualCorrections.png'];
        saveas(h, currentFile);
        currentFrame = imread(currentFile);
        imwrite(currentFrame, currentFile, 'png');
        
        save([configRoot filesep 'driftDatabase.mat'], 'smoothedDriftError', '-append');
    end
    
    figure;
    plot(timepoints, finalDriftOffsets);
    title('Final offset correction');
    legend('xOffset', 'yOffset', 'zOffset');
    xlabel('Time point');
    ylabel('Offset (pixel)');
    h = gca;
    currentFile = [configRoot filesep 'driftStep6_finalOffsets.png'];
    saveas(h, currentFile);
    currentFrame = imread(currentFile);
    imwrite(currentFrame, currentFile, 'png');
end;