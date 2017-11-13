%%% Generation of single-view weight matrices for MVD %%%

% Conceptual processing workflow:

% 1) Read specimen tracking text file and convert to Matlab matrix
%    (Time, Specimen, Camera, Type, Backgrd, XYZ_Stage_Drift, #newPlnLow, #newPlnHi)

% 2) Reconstruct artificial (x, y, z)-drift vector as a function of time by combining
%    specimen tracking corrections (XYZ_Stage_Drift) and z-range expansion corrections
%    (#newPlnLow, #newPlnHi), signifying total artificial voxel shift measured relative
%    to current stack origin (0, 0, 0)

% 3) Read image stacks (four views: SPM00_CM00_CHN00, SPM00_CM01_CHN00, SPM01_CM00_CHN00,
%    SPM01_CM01_CHN00), apply adaptive threshold (Gaussian smoothing, threshold proportional
%    to mean intensity; following approach used in clusterPT/clusterMF) and calculate
%    centroids of binary stacks

% 4) Determine average centroids for each specimen (compute average of CM0 and CM1 centroids
%    in order to obtain estimate of specimen centroid, which represents the coordinates
%    relevant for masking)

% 5) Subtract artificial drift vectors from average centroids determined in step (4)
%    and compute smoothed version of these baseline-corrected centroid trajectories
%    (i.e. smoothing is performed in sample space)

% 6) Add artificial drift vectors back to smoothed version of baseline-corrected
%    centroid trajectories (i.e. smoothed data is converted back to image space)

% 7) Compute masks with �fading zones� (smooth transition from 1 to 0, possibly
%    testing linear/sigmoidal/exponential transitions) positioned relative to centroids
%    determined in step (6), thus producing masks that allow us to prevent the contribution
%    of blurry regions in each view to MVD results

%% Configuration

dataFolder       = ['S:' filesep 'SiMView1' filesep '15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected']; % image data input folder
metaFolder       = ['S:' filesep 'SiMView1' filesep '15-08-10' filesep 'SpecimenTrackingFiles' filesep 'Tracking_160708'];   % meta data input folder (SiMView specimen tracking)
                                                                                   % set to '' if specimen tracking and z-range exansion logs are not available
inputIdentifier  = 'unmasked_r5_d3_e2';                                            % identification string used in file names for input data
outputIdentifier = 'unmasked_r5_d3_e2.threshold_0.2';                              % identification string used in file names for output data

voxelSize        = [6.5/16 2.031];

specimen         = 0;
timepoints       = 0:1243;
channel          = 0;

restrictXRange   = [0 1 2048]; % flag for restricting x-range of image stack used for mask computation (first value)
                               % if first value is set to 1, the second and third values indicate start and end x-coordinates of the ROI
splitXRange      = 1;          % construct transition zones parallel to the x-y plane (value 0) or tilted around the y-axis (value 1)
splitting        = 5;          % splitting factor for reducing memory consumption for imgaussianAnisotropy
gaussianKernel   = [7 3];      % size and sigma for Gaussian smoothing
checkForZeros    = 0;          % only consider non-zero voxels for background computation (set to zero if stack is unmasked)
percentile       = [1 100];    % percentile and subsampling for background computation
threshold        = 0.2;        % adaptive threshold for foreground computation

smoothing        = [0 8];      % smoothing flag ('rloess'), smoothing window size
smoothingStart   = '';         % only relevant if smoothing(1) == 1:
                               % start time point index for smoothing operation, set to '' if entire timeline should be smoothed
subtractBaseline = 1;          % only relevant if smoothing(1) == 1:
                               % flag for decoupling data from drift correction and z-range expansion before smoothing centroid trajectories,
                               % baseline will be added back to the trajectories after smoothing
lockLowCentroid  = [0 1050];   % flag for keeping lower centroid at a fixed distance to upper centroid (use last distance), start time index for fixing distance
                               % note: currently only implemented for isempty(metaFolder)==1, but can in principle be generally applied
primaryCamera    = 0;          % 0: high-contrast information is stored in low-z planes of first camera (camera 0)
                               % 1: high-contrast information is stored in low-z planes of second camera (camera 1)
transitionDelta  = 10;         % distance of (start point of) transition zone from centroid z-location, unit: z planes
transitionLength = 20;         % size of transition zone (distance across which mask values drop from 1 to 0), unit: z planes
transitionType   = [1 1.1];    % 0: linear transition,    (start+transitionLength-z)/transitionLength with z = start .. start+transitionLength
                               % 1: sigmoidal transition, (1/(1+exp(-(z-center)*10/(width*transitionLength))) with z = center-transitionLength/2 .. center+transitionLength/2
                               %    width modifier is provided as a second value
transitionToZero = [0 100];    % 0: fade to background level defined during masking (if second value is 0) or to second value (if second value is >0)
                               % 1: fade to absolute zero intensity

forceCentroids   = 0;          % flag for enforcing centroid computations even if existing results are detected
forceWeights     = 0;          % flag for enforcing weights matrix computations even if existing results are detected

saveMasks        = 1;          % flag for saving binary masks used to compute centroids
applyWeights     = 1;          % flag for computing and saving weighted image stacks
saveProjections  = 1;          % only relevant if applyWeights == 1:
                               % flag for computing and saving xy/xz/yz projections of weighted image stacks

poolWorkers      = 12;         % number of parallel workers

%% Analyze experiment meta data

if ~isempty(inputIdentifier)
    inputIdentifier = ['.' inputIdentifier];
end;

if ~isempty(outputIdentifier)
    outputIdentifier = ['_' outputIdentifier];
end;

if ~isempty(metaFolder)
    disp(' ');
    
    metaDataFilename = [metaFolder filesep 'MetaMatrix.mat'];
    
    if exist(metaDataFilename, 'file')~=2
        disp('Parsing specimen tracking log file');
        
        trackingLogFilename = [metaFolder filesep 'Specimen tracking diagnostics.txt'];
        
        fid = fopen(trackingLogFilename, 'r');
        if fid<1
            error('ERROR: Could not open file %s for reading', filename);
        end;
        
        metaData = [];
        while true
            s = fgetl(fid);
            if ~ischar(s)
                break;
            end;
            m = regexp(s, '\t', 'split');
            if numel(m) == 15 && ~strcmp(m{1}, 'Time') && strcmp(m{3}, '0') && strcmp(m{4}, 'X/Z')
                currentTime = str2double(m{1});
                currentSpecimen = str2double(m{2});
                currentNewPlanesLow = str2double(m{11});
                currentNewPlanesHigh = str2double(m{12});
                
                driftVector = regexp(m{5}, '_', 'split');
                xComponents = regexp(driftVector{1}, '=', 'split');
                yComponents = regexp(driftVector{2}, '=', 'split');
                zComponents = regexp(driftVector{3}, '=', 'split');
                xDrift = str2double(xComponents{2});
                yDrift = str2double(yComponents{2});
                zDrift = str2double(zComponents{2});
                
                metaData = cat(1, metaData, [currentTime currentSpecimen, xDrift, yDrift, zDrift, currentNewPlanesLow, currentNewPlanesHigh]);
            end;
        end;
        fclose(fid);
        
        disp('Validating z-range expansion meta data');
        
        discrepanciesDetected = 0;
        for i = 1:size(metaData, 1)
            if (metaData(i, 6)~=0 || metaData(i, 7)~=0) && i<(size(metaData, 1)-1)
                currentT = metaData(i, 1);
                currentS = metaData(i, 2);
                stackName1Camera0 = [dataFolder filesep 'SPM' num2str(currentS, '%.2d') filesep 'TM' num2str(currentT, '%.6d') filesep 'SPM' num2str(currentS, '%.2d') ...
                    '_TM' num2str(currentT, '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
                stackName2Camera0 = [dataFolder filesep 'SPM' num2str(currentS, '%.2d') filesep 'TM' num2str(currentT+1, '%.6d') filesep 'SPM' num2str(currentS, '%.2d') ...
                    '_TM' num2str(currentT+1, '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
                stack1Camera0Header = readKLBheader(stackName1Camera0);
                stack2Camera0Header = readKLBheader(stackName2Camera0);
                camera0zDifference = stack2Camera0Header.xyzct(3) - stack1Camera0Header.xyzct(3);
                if camera0zDifference ~= (metaData(i, 6)+metaData(i, 7))
                    warning(['z-range expansion meta data does not match image data for camera 0 of specimen ' num2str(currentS) ' at time point ' num2str(currentT)]);
                    discrepanciesDetected = discrepanciesDetected+1;
                end;
                stackName1Camera1 = [dataFolder filesep 'SPM' num2str(currentS, '%.2d') filesep 'TM' num2str(currentT, '%.6d') filesep 'SPM' num2str(currentS, '%.2d') ...
                    '_TM' num2str(currentT, '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
                stackName2Camera1 = [dataFolder filesep 'SPM' num2str(currentS, '%.2d') filesep 'TM' num2str(currentT+1, '%.6d') filesep 'SPM' num2str(currentS, '%.2d') ...
                    '_TM' num2str(currentT+1, '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
                stack1Camera1Header = readKLBheader(stackName1Camera1);
                stack2Camera1Header = readKLBheader(stackName2Camera1);
                camera1zDifference = stack2Camera1Header.xyzct(3) - stack1Camera1Header.xyzct(3);
                if camera1zDifference ~= (metaData(i, 6)+metaData(i, 7))
                    warning(['z-range expansion meta data does not match image data for camera 1 of specimen ' num2str(currentS) ' at time point ' num2str(currentT)]);
                    discrepanciesDetected = discrepanciesDetected+1;
                end;
            end;
        end;
        specimen0AddedPlanes = sum(metaData(metaData(1:(end-2), 2) == 0, 6)) + sum(metaData(metaData(1:(end-2), 2) == 0, 7));
        specimen1AddedPlanes = sum(metaData(metaData(1:(end-2), 2) == 1, 6)) + sum(metaData(metaData(1:(end-2), 2) == 1, 7));
        stackName1Specimen0Camera0 = [dataFolder filesep 'SPM00' filesep 'TM' num2str(min(metaData(:, 1)), '%.6d') filesep 'SPM00_TM' num2str(min(metaData(:, 1)), '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stackName2Specimen0Camera0 = [dataFolder filesep 'SPM00' filesep 'TM' num2str(max(metaData(:, 1)), '%.6d') filesep 'SPM00_TM' num2str(max(metaData(:, 1)), '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stack1Specimen0Camera0Header = readKLBheader(stackName1Specimen0Camera0);
        stack2Specimen0Camera0Header = readKLBheader(stackName2Specimen0Camera0);
        specimen0Camera0zDifference = stack2Specimen0Camera0Header.xyzct(3) - stack1Specimen0Camera0Header.xyzct(3);
        stackName1Specimen0Camera1 = [dataFolder filesep 'SPM00' filesep 'TM' num2str(min(metaData(:, 1)), '%.6d') filesep 'SPM00_TM' num2str(min(metaData(:, 1)), '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stackName2Specimen0Camera1 = [dataFolder filesep 'SPM00' filesep 'TM' num2str(max(metaData(:, 1)), '%.6d') filesep 'SPM00_TM' num2str(max(metaData(:, 1)), '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stack1Specimen0Camera1Header = readKLBheader(stackName1Specimen0Camera1);
        stack2Specimen0Camera1Header = readKLBheader(stackName2Specimen0Camera1);
        specimen0Camera1zDifference = stack2Specimen0Camera1Header.xyzct(3) - stack1Specimen0Camera1Header.xyzct(3);
        stackName1Specimen1Camera0 = [dataFolder filesep 'SPM01' filesep 'TM' num2str(min(metaData(:, 1)), '%.6d') filesep 'SPM01_TM' num2str(min(metaData(:, 1)), '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stackName2Specimen1Camera0 = [dataFolder filesep 'SPM01' filesep 'TM' num2str(max(metaData(:, 1)), '%.6d') filesep 'SPM01_TM' num2str(max(metaData(:, 1)), '%.6d') '_CM00_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stack1Specimen1Camera0Header = readKLBheader(stackName1Specimen1Camera0);
        stack2Specimen1Camera0Header = readKLBheader(stackName2Specimen1Camera0);
        specimen1Camera0zDifference = stack2Specimen1Camera0Header.xyzct(3) - stack1Specimen1Camera0Header.xyzct(3);
        stackName1Specimen1Camera1 = [dataFolder filesep 'SPM01' filesep 'TM' num2str(min(metaData(:, 1)), '%.6d') filesep 'SPM01_TM' num2str(min(metaData(:, 1)), '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stackName2Specimen1Camera1 = [dataFolder filesep 'SPM01' filesep 'TM' num2str(max(metaData(:, 1)), '%.6d') filesep 'SPM01_TM' num2str(max(metaData(:, 1)), '%.6d') '_CM01_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
        stack1Specimen1Camera1Header = readKLBheader(stackName1Specimen1Camera1);
        stack2Specimen1Camera1Header = readKLBheader(stackName2Specimen1Camera1);
        specimen1Camera1zDifference = stack2Specimen1Camera1Header.xyzct(3) - stack1Specimen1Camera1Header.xyzct(3);
        if (specimen0AddedPlanes~=specimen0Camera0zDifference) || (specimen0AddedPlanes~=specimen0Camera1zDifference)
            warning('z-range expansion meta data for specimen 0 does not agree with image data across the full time range');
            discrepanciesDetected = discrepanciesDetected+1;
        end;
        if (specimen1AddedPlanes~=specimen1Camera0zDifference) || (specimen1AddedPlanes~=specimen1Camera1zDifference)
            warning('z-range expansion meta data for specimen 1 does not agree with image data across the full time range');
            discrepanciesDetected = discrepanciesDetected+1;
        end;
        
        if discrepanciesDetected>0
            error(['ERROR: ' num2str(discrepanciesDetected) ' instances of a discrepancy between meta data and image data detected - aborting processing']);
        else
            save(metaDataFilename, 'metaData');
            disp('Meta data matrix passed consistency checks and was saved to disk');
        end;
    else
        load(metaDataFilename);
        disp('Existing meta data matrix found and loaded');
    end;
    
    metaData = metaData(metaData(:, 2)==specimen, [1 3:7]);
    metaData(:, 2) = metaData(:, 2)/voxelSize(1);
    metaData(:, 3) = metaData(:, 3)/voxelSize(1);
    metaData(:, 4) = metaData(:, 4)/voxelSize(2);
    processedMetaData = cat(2, metaData(:, 1:4), zeros(size(metaData, 1), 4));
    for t = 1:(size(processedMetaData, 1)-1) % only NewLowPlanes are relevant for measurements relative to z=0
        processedMetaData((t+1):end, 5) = processedMetaData((t+1):end, 5) + metaData(t, 5); % cumulative shift resulting from NewPlanesLow
        processedMetaData((t+1):end, 6) = processedMetaData((t+1):end, 6) + metaData(t, 6); % cumulative shift resulting from NewPlanesHigh
        processedMetaData((t+1):end, 7) = processedMetaData((t+1):end, 7) + metaData(t, 5); % combined effects of drift and z-range expansion (additive)
        processedMetaData((t+1):end, 8) = processedMetaData((t+1):end, 8) - metaData(t, 5); % combined effects of drift and z-range expansion (subtractive)
    end;
    
    xDriftCorrections = [];
    yDriftCorrections = [];
    zDriftCorrections = [];
    zLowPlaneAdditions = [];
    zHighPlaneAdditions = [];
    for t = 1:numel(timepoints)
        xDriftCorrections = cat(1, xDriftCorrections, processedMetaData(processedMetaData(:, 1) == timepoints(t), 2));
        yDriftCorrections = cat(1, yDriftCorrections, processedMetaData(processedMetaData(:, 1) == timepoints(t), 3));
        zDriftCorrections = cat(1, zDriftCorrections, processedMetaData(processedMetaData(:, 1) == timepoints(t), 4));
        zLowPlaneAdditions = cat(1, zLowPlaneAdditions, processedMetaData(processedMetaData(:, 1) == timepoints(t), 5));
        zHighPlaneAdditions = cat(1, zHighPlaneAdditions, processedMetaData(processedMetaData(:, 1) == timepoints(t), 6));
    end;
end;

%% Determine foreground centroids

kernelSize  = [gaussianKernel(1) gaussianKernel(1) max(1, gaussianKernel(1)/(voxelSize(2)/voxelSize(1)))];
kernelSigma = [gaussianKernel(2) gaussianKernel(2) max(1, gaussianKernel(2)/(voxelSize(2)/voxelSize(1)))];

xSizes             = zeros(numel(timepoints), 2);
ySizes             = zeros(numel(timepoints), 2);
zSizes             = zeros(numel(timepoints), 2);
minIntensities     = zeros(numel(timepoints), 2);
meanIntensities    = zeros(numel(timepoints), 2);
absoluteThresholds = zeros(numel(timepoints), 2);
if splitXRange
    centroidsX     = zeros(numel(timepoints), 2, 2);
    centroidsY     = zeros(numel(timepoints), 2, 2);
    centroidsZ     = zeros(numel(timepoints), 2, 2);
else
    centroidsX     = zeros(numel(timepoints), 2);
    centroidsY     = zeros(numel(timepoints), 2);
    centroidsZ     = zeros(numel(timepoints), 2);
end;

disp(' ');
disp('Determining foreground centroids');

for t = 1:numel(timepoints)
    disp(['*** Processing time point ' num2str(timepoints(t))]);
    
    for camera = [0 1]
        if restrictXRange(1) && splitXRange
            xRangeString = ['_X' num2str(restrictXRange(2)) '-' num2str(restrictXRange(3)) '-Split'];
        elseif restrictXRange(1) && ~splitXRange
            xRangeString = ['_X' num2str(restrictXRange(2)) '-' num2str(restrictXRange(3))];
        elseif ~restrictXRange(1) && splitXRange
            xRangeString = '_XSplit';
        else
            xRangeString = '';
        end;
        maskFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') inputIdentifier '.binaryMask_' num2str(threshold) xRangeString '.klb'];
        logsFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') inputIdentifier '.processingLogs_' num2str(threshold) xRangeString '.mat'];
        
        if exist(logsFilename, 'file')~=2 || forceCentroids
            disp(['    Creating binary mask for camera ' num2str(camera)]);
            
            stackFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
                filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
            stack = readImage(stackFilename);
            
            [xSize, ySize, zSize] = size(stack);
            xSizes(t, camera+1) = xSize;
            ySizes(t, camera+1) = ySize;
            zSizes(t, camera+1) = zSize;
            
            if restrictXRange(1)
                stack = stack(restrictXRange(2):restrictXRange(3), :, :);
                xSizeForProcessing = size(stack, 1);
            else
                xSizeForProcessing = xSize;
            end;
            
            if checkForZeros
                nonZeroStack = stack(stack>0);
                minIntensity = prctile(nonZeroStack(1:percentile(2):end), percentile(1));
                clear nonZeroStack;
            else
                minIntensity = prctile(stack(1:percentile(2):end), percentile(1));
            end;
            minIntensities(t, camera+1) = minIntensity;
            
            if splitting>1
                gaussStack = zeros(xSizeForProcessing, ySize, zSize, 'uint16');
                splittingMargin = 2*gaussianKernel(1);
                
                for i = 1:splitting
                    slabStart = max(1, round(((i-1)*xSizeForProcessing/splitting)+1-splittingMargin));
                    slabStop = min(xSizeForProcessing, round((i*xSizeForProcessing/splitting)+splittingMargin));
                    convolvedSlab = uint16(imgaussianAnisotropy(double(stack(slabStart:slabStop, :, :)), kernelSigma, kernelSize));
                    if i==1
                        gaussStack(1:(slabStop-splittingMargin), :, :) = convolvedSlab(1:(end-splittingMargin), :, :);
                    elseif i==splitting
                        gaussStack((slabStart+splittingMargin):end, :, :) = convolvedSlab((1+splittingMargin):end, :, :);
                    else % i > 1 && i < splitting
                        gaussStack((slabStart+splittingMargin):(slabStop-splittingMargin), :, :) = convolvedSlab((1+splittingMargin):(end-splittingMargin), :, :);
                    end;
                    clear convolvedSlab;
                end;
            else
                gaussStack = uint16(imgaussianAnisotropy(double(stack), kernelSigmaArray, kernelSizeArray));
            end;
            clear stack;
            
            intensityStatistics = zeros(splitting, 2);
            for i = 1:splitting
                slabStart = round(((i-1)*xSizeForProcessing/splitting)+1);
                slabStop  = round(i*xSizeForProcessing/splitting);
                temporaryArray = gaussStack(slabStart:slabStop, :, :);
                temporaryArray = temporaryArray(temporaryArray > minIntensity);
                intensityStatistics(i, 1) = sum(temporaryArray(:));
                intensityStatistics(i, 2) = size(temporaryArray, 1);
            end;
            meanIntensity = sum(intensityStatistics(:, 1))/sum(intensityStatistics(:, 2));
            clear temporaryArray intensityStatistics;
            meanIntensities(t, camera+1) = meanIntensity;
            
            absoluteThreshold = minIntensity+(meanIntensity-minIntensity)*threshold;
            absoluteThresholds(t, camera+1) = absoluteThreshold;
            
            mask = gaussStack>absoluteThreshold;
            if saveMasks
                writeImage(uint8(mask)*255, maskFilename);
            end;
            clear gaussStack;
            
            disp(['    Determining centroid coordinates for camera ' num2str(camera)]);
            
            if splitXRange
                xCenter = round(xSizeForProcessing/2);
                centroidX = zeros(1, 2);
                centroidY = zeros(1, 2);
                centroidZ = zeros(1, 2);
                
                coordinateMask = repmat((1:xCenter)', [1 ySize zSize]);
                coordinateMask(mask(1:xCenter,:,:)==0) = NaN;
                centroidX(1) = nanmean(coordinateMask(:));
                centroidsX(t, camera+1, 1) = centroidX(1);
                coordinateMask = repmat(((xCenter+1):xSizeForProcessing)', [1 ySize zSize]);
                coordinateMask(mask((xCenter+1):xSizeForProcessing,:,:)==0) = NaN;
                centroidX(2) = nanmean(coordinateMask(:));
                centroidsX(t, camera+1, 2) = centroidX(2);
                
                coordinateMask = repmat(1:ySize, [xCenter 1 zSize]);
                coordinateMask(mask(1:xCenter,:,:)==0) = NaN;
                centroidY(1) = nanmean(coordinateMask(:));
                centroidsY(t, camera+1, 1) = centroidY(1);
                coordinateMask = repmat(1:ySize, [(xSizeForProcessing-xCenter) 1 zSize]);
                coordinateMask(mask((xCenter+1):xSizeForProcessing,:,:)==0) = NaN;
                centroidY(2) = nanmean(coordinateMask(:));
                centroidsY(t, camera+1, 2) = centroidY(2);
                
                coordinateMask = repmat(reshape(1:zSize, [1 1 zSize]), [xCenter ySize 1]);
                coordinateMask(mask(1:xCenter,:,:)==0) = NaN;
                centroidZ(1) = nanmean(coordinateMask(:));
                centroidsZ(t, camera+1, 1) = centroidZ(1);
                coordinateMask = repmat(reshape(1:zSize, [1 1 zSize]), [(xSizeForProcessing-xCenter) ySize 1]);
                coordinateMask(mask((xCenter+1):xSizeForProcessing,:,:)==0) = NaN;
                centroidZ(2) = nanmean(coordinateMask(:));
                centroidsZ(t, camera+1, 2) = centroidZ(2);
            else
                coordinateMask = repmat((1:xSizeForProcessing)', [1 ySize zSize]);
                coordinateMask(mask == 0) = NaN;
                centroidX = nanmean(coordinateMask(:));
                centroidsX(t, camera+1) = centroidX;
                
                coordinateMask = repmat(1:ySize, [xSizeForProcessing 1 zSize]);
                coordinateMask(mask == 0) = NaN;
                centroidY = nanmean(coordinateMask(:));
                centroidsY(t, camera+1) = centroidY;
                
                coordinateMask = repmat(reshape(1:zSize, [1 1 zSize]), [xSizeForProcessing ySize 1]);
                coordinateMask(mask == 0) = NaN;
                centroidZ = nanmean(coordinateMask(:));
                centroidsZ(t, camera+1) = centroidZ;
            end;
            clear coordinateMask;
            
            save(logsFilename, 'xSize', 'ySize', 'zSize', 'minIntensity', 'meanIntensity', 'absoluteThreshold', 'centroidX', 'centroidY', 'centroidZ');
        else
            disp(['    Loading processing logs for camera ' num2str(camera)]);
            
            load(logsFilename);

            xSizes(t, camera+1) = xSize;
            ySizes(t, camera+1) = ySize;
            zSizes(t, camera+1) = zSize;
            minIntensities(t, camera+1) = minIntensity;
            meanIntensities(t, camera+1) = meanIntensity;
            absoluteThresholds(t, camera+1) = absoluteThreshold;
            if splitXRange
                centroidsX(t, camera+1, 1) = centroidX(1);
                centroidsY(t, camera+1, 1) = centroidY(1);
                centroidsZ(t, camera+1, 1) = centroidZ(1);
                centroidsX(t, camera+1, 2) = centroidX(2);
                centroidsY(t, camera+1, 2) = centroidY(2);
                centroidsZ(t, camera+1, 2) = centroidZ(2);
            else
                centroidsX(t, camera+1) = centroidX;
                centroidsY(t, camera+1) = centroidY;
                centroidsZ(t, camera+1) = centroidZ;
            end;
        end;
    end;
end;

%% Process centroid trajectories

disp(' ');
disp('Processing centroid trajectories');

if ~isempty(metaFolder)
    if splitXRange
        averagedCentroidsZ1 = 0.5*centroidsZ(:, 1, 1)+0.5*centroidsZ(:, 2, 1);
        averagedCentroidsZ2 = 0.5*centroidsZ(:, 1, 2)+0.5*centroidsZ(:, 2, 2);
        
        centroidsRawPlotsFileName1 = ['createMasks.centroidsRaw_SPM' num2str(specimen) outputIdentifier '_X1.png'];
        centroidsRawPlotsFileName2 = ['createMasks.centroidsRaw_SPM' num2str(specimen) outputIdentifier '_X2.png'];
        centroidsBaselineCorrectedPlotsFileName1 = ['createMasks.centroidsBaselineCorrected_SPM' num2str(specimen) outputIdentifier '_X1.png'];
        centroidsBaselineCorrectedPlotsFileName2 = ['createMasks.centroidsBaselineCorrected_SPM' num2str(specimen) outputIdentifier '_X2.png'];
        centroidsSmoothedPlotsFileName1 = ['createMasks.centroidsSmoothed_SPM' num2str(specimen) outputIdentifier '_X1.png'];
        centroidsSmoothedPlotsFileName2 = ['createMasks.centroidsSmoothed_SPM' num2str(specimen) outputIdentifier '_X2.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ1, 'k', 'LineWidth', 1);
        hold on;
        plot(timepoints, zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1);
        plot(timepoints, zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        title('Raw image centroid, drift correction and z-range expansion data');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids (top half)', 'Drift correction', 'Z-range expansion (low z only)', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsRawPlotsFileName1);
        
        figure;
        plot(timepoints, averagedCentroidsZ2, 'k', 'LineWidth', 1);
        hold on;
        plot(timepoints, zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1);
        plot(timepoints, zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        title('Raw image centroid, drift correction and z-range expansion data');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids (bottom half)', 'Drift correction', 'Z-range expansion (low z only)', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsRawPlotsFileName2);
        
        figure;
        plot(timepoints, averagedCentroidsZ1, 'k', 'LineWidth', 1); % Effect 1: z-drift correction in +z (stage interval shifted to higher z)
        %           shifts centroids artificially to lower z and vice versa
        %           -> add z-drift vector to obtain real [smooth] centroid motion
        % Effect 2: adding low-z planes shifts centroids artifically to higher z
        %           -> subtract low-z planes to obtain real [smooth] centroid motion
        % Combined baseline correction: (1) add "zDriftCorrections", (2) subtract "zLowPlaneAdditions"
        hold on;
        plot(timepoints, averagedCentroidsZ1+zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1)
        plot(timepoints, averagedCentroidsZ1-zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        plot(timepoints, averagedCentroidsZ1+zDriftCorrections-zLowPlaneAdditions, 'Color', [0 0 1], 'LineWidth', 1)
        title('Baseline subtraction of centroid trajectories');
        xlabel('Time point');
        ylabel('Centroid z-position (top half, unit: z-planes)');
        legend('Raw', 'Drift decoupled', 'Expansion decoupled', 'Drift and expansion decoupled', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsBaselineCorrectedPlotsFileName1);
        
        figure;
        plot(timepoints, averagedCentroidsZ2, 'k', 'LineWidth', 1); % Effect 1: z-drift correction in +z (stage interval shifted to higher z)
        %           shifts centroids artificially to lower z and vice versa
        %           -> add z-drift vector to obtain real [smooth] centroid motion
        % Effect 2: adding low-z planes shifts centroids artifically to higher z
        %           -> subtract low-z planes to obtain real [smooth] centroid motion
        % Combined baseline correction: (1) add "zDriftCorrections", (2) subtract "zLowPlaneAdditions"
        hold on;
        plot(timepoints, averagedCentroidsZ2+zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1)
        plot(timepoints, averagedCentroidsZ2-zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        plot(timepoints, averagedCentroidsZ2+zDriftCorrections-zLowPlaneAdditions, 'Color', [0 0 1], 'LineWidth', 1)
        title('Baseline subtraction of centroid trajectories');
        xlabel('Time point');
        ylabel('Centroid z-position (bottom half, unit: z-planes)');
        legend('Raw', 'Drift decoupled', 'Expansion decoupled', 'Drift and expansion decoupled', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsBaselineCorrectedPlotsFileName2);
        
        if subtractBaseline && smoothing(1)
            correctedCentroidsZ1 = averagedCentroidsZ1+zDriftCorrections-zLowPlaneAdditions;
            
            if ~isempty(smoothingStart)
                smoothCentroidsZ1 = correctedCentroidsZ1;
                smoothCentroidsZ1(smoothingStart:end) = smooth(correctedCentroidsZ1(smoothingStart:end), smoothing(2), 'rloess');
            else
                smoothCentroidsZ1 = smooth(correctedCentroidsZ1, smoothing(2), 'rloess');
            end;
            
            finalCentroidsZ1 = smoothCentroidsZ1-zDriftCorrections+zLowPlaneAdditions;
            processedCentroidsZ1 = [averagedCentroidsZ1 correctedCentroidsZ1 smoothCentroidsZ1 finalCentroidsZ1];
            
            figure;
            plot(timepoints, averagedCentroidsZ1, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, correctedCentroidsZ1, 'Color', [1 0.5 0.5], 'LineWidth', 1);
            plot(timepoints, smoothCentroidsZ1, 'Color', [0.5 0 0], 'LineWidth', 1);
            plot(timepoints, finalCentroidsZ1, 'Color', [0 0 0.5], 'LineWidth', 1);
        elseif smoothing(1)
            if ~isempty(smoothingStart)
                finalCentroidsZ1 = averagedCentroidsZ1;
                finalCentroidsZ1(smoothingStart:end) = smooth(averagedCentroidsZ1(smoothingStart:end), smoothing(2), 'rloess');
            else
                finalCentroidsZ1 = smooth(averagedCentroidsZ1, smoothing(2), 'rloess');
            end;
            
            processedCentroidsZ1 = [averagedCentroidsZ1 finalCentroidsZ1];
            
            figure;
            plot(timepoints, averagedCentroidsZ1, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, finalCentroidsZ1, 'Color', [0 0 0.5], 'LineWidth', 1);
        else
            finalCentroidsZ1 = averagedCentroidsZ1;
            processedCentroidsZ1 = finalCentroidsZ1;
        end;
        
        if smoothing(1)
            title('Smooth centroid trajectories');
            xlabel('Time point');
            ylabel('Centroid z-position (top half, unit: z-planes)');
            if subtractBaseline
                legend('Raw', 'Baseline-subtracted', 'Baseline-subtracted smooth', 'Baseline-reconstituted smooth (final)', 'Location', 'SouthOutside');
            else
                legend('Raw', 'Smooth (final)', 'Location', 'SouthOutside');
            end;
            h = gca;
            saveas(h, centroidsSmoothedPlotsFileName1);
        end;
        
        if subtractBaseline && smoothing(1)
            correctedCentroidsZ2 = averagedCentroidsZ2+zDriftCorrections-zLowPlaneAdditions;
            
            if ~isempty(smoothingStart)
                smoothCentroidsZ2 = correctedCentroidsZ2;
                smoothCentroidsZ2(smoothingStart:end) = smooth(correctedCentroidsZ2(smoothingStart:end), smoothing(2), 'rloess');
            else
                smoothCentroidsZ2 = smooth(correctedCentroidsZ2, smoothing(2), 'rloess');
            end;
            
            finalCentroidsZ2 = smoothCentroidsZ2-zDriftCorrections+zLowPlaneAdditions;
            processedCentroidsZ2 = [averagedCentroidsZ2 correctedCentroidsZ2 smoothCentroidsZ2 finalCentroidsZ2];
            
            figure;
            plot(timepoints, averagedCentroidsZ2, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, correctedCentroidsZ2, 'Color', [1 0.5 0.5], 'LineWidth', 1);
            plot(timepoints, smoothCentroidsZ2, 'Color', [0.5 0 0], 'LineWidth', 1);
            plot(timepoints, finalCentroidsZ2, 'Color', [0 0 0.5], 'LineWidth', 1);
        elseif smoothing(1)
            if ~isempty(smoothingStart)
                finalCentroidsZ2 = averagedCentroidsZ2;
                finalCentroidsZ2(smoothingStart:end) = smooth(averagedCentroidsZ2(smoothingStart:end), smoothing(2), 'rloess');
            else
                finalCentroidsZ2 = smooth(averagedCentroidsZ2, smoothing(2), 'rloess');
            end;
            
            processedCentroidsZ2 = [averagedCentroidsZ2 finalCentroidsZ2];
            
            figure;
            plot(timepoints, averagedCentroidsZ2, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, finalCentroidsZ2, 'Color', [0 0 0.5], 'LineWidth', 1);
        else
            finalCentroidsZ2 = averagedCentroidsZ2;
            processedCentroidsZ2 = finalCentroidsZ2;
        end;
        
        if smoothing(1)
            title('Smooth centroid trajectories');
            xlabel('Time point');
            ylabel('Centroid z-position (bottom half, unit: z-planes)');
            if subtractBaseline
                legend('Raw', 'Baseline-subtracted', 'Baseline-subtracted smooth', 'Baseline-reconstituted smooth (final)', 'Location', 'SouthOutside');
            else
                legend('Raw', 'Smooth (final)', 'Location', 'SouthOutside');
            end;
            h = gca;
            saveas(h, centroidsSmoothedPlotsFileName2);
        end;
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'processedMetaData', 'xDriftCorrections', 'yDriftCorrections', 'zDriftCorrections', 'zLowPlaneAdditions', 'zHighPlaneAdditions', ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ1', 'processedCentroidsZ2');
    else
        averagedCentroidsZ = 0.5*centroidsZ(:, 1)+0.5*centroidsZ(:, 2);
        
        centroidsRawPlotsFileName = ['createMasks.centroidsRaw_SPM' num2str(specimen) outputIdentifier '.png'];
        centroidsBaselineCorrectedPlotsFileName = ['createMasks.centroidsBaselineCorrected_SPM' num2str(specimen) outputIdentifier '.png'];
        centroidsSmoothedPlotsFileName = ['createMasks.centroidsSmoothed_SPM' num2str(specimen) outputIdentifier '.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ, 'k', 'LineWidth', 1);
        hold on;
        plot(timepoints, zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1);
        plot(timepoints, zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        title('Raw image centroid, drift correction and z-range expansion data');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids', 'Drift correction', 'Z-range expansion (low z only)', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsRawPlotsFileName);
        
        figure;
        plot(timepoints, averagedCentroidsZ, 'k', 'LineWidth', 1); % Effect 1: z-drift correction in +z (stage interval shifted to higher z)
        %           shifts centroids artificially to lower z and vice versa
        %           -> add z-drift vector to obtain real [smooth] centroid motion
        % Effect 2: adding low-z planes shifts centroids artifically to higher z
        %           -> subtract low-z planes to obtain real [smooth] centroid motion
        % Combined baseline correction: (1) add "zDriftCorrections", (2) subtract "zLowPlaneAdditions"
        hold on;
        plot(timepoints, averagedCentroidsZ+zDriftCorrections, 'Color', [0.5 0 0], 'LineWidth', 1)
        plot(timepoints, averagedCentroidsZ-zLowPlaneAdditions, 'Color', [0 0.5 0], 'LineWidth', 1);
        plot(timepoints, averagedCentroidsZ+zDriftCorrections-zLowPlaneAdditions, 'Color', [0 0 1], 'LineWidth', 1)
        title('Baseline subtraction of centroid trajectories');
        xlabel('Time point');
        ylabel('Centroid z-position (unit: z-planes)');
        legend('Raw', 'Drift decoupled', 'Expansion decoupled', 'Drift and expansion decoupled', 'Location', 'SouthOutside');
        h = gca;
        saveas(h, centroidsBaselineCorrectedPlotsFileName);
        
        if subtractBaseline && smoothing(1)
            correctedCentroidsZ = averagedCentroidsZ+zDriftCorrections-zLowPlaneAdditions;
            
            if ~isempty(smoothingStart)
                smoothCentroidsZ = correctedCentroidsZ;
                smoothCentroidsZ(smoothingStart:end) = smooth(correctedCentroidsZ(smoothingStart:end), smoothing(2), 'rloess');
            else
                smoothCentroidsZ = smooth(correctedCentroidsZ, smoothing(2), 'rloess');
            end;
            
            finalCentroidsZ = smoothCentroidsZ-zDriftCorrections+zLowPlaneAdditions;
            processedCentroidsZ = [averagedCentroidsZ correctedCentroidsZ smoothCentroidsZ finalCentroidsZ];
            
            figure;
            plot(timepoints, averagedCentroidsZ, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, correctedCentroidsZ, 'Color', [1 0.5 0.5], 'LineWidth', 1);
            plot(timepoints, smoothCentroidsZ, 'Color', [0.5 0 0], 'LineWidth', 1);
            plot(timepoints, finalCentroidsZ, 'Color', [0 0 0.5], 'LineWidth', 1);
        elseif smoothing(1)
            if ~isempty(smoothingStart)
                finalCentroidsZ = averagedCentroidsZ;
                finalCentroidsZ(smoothingStart:end) = smooth(averagedCentroidsZ(smoothingStart:end), smoothing(2), 'rloess');
            else
                finalCentroidsZ = smooth(averagedCentroidsZ, smoothing(2), 'rloess');
            end;
            
            processedCentroidsZ = [averagedCentroidsZ finalCentroidsZ];
            
            figure;
            plot(timepoints, averagedCentroidsZ, 'Color', [0.5 0.5 1], 'LineWidth', 1);
            hold on;
            plot(timepoints, finalCentroidsZ, 'Color', [0 0 0.5], 'LineWidth', 1);
        else
            finalCentroidsZ = averagedCentroidsZ;
            processedCentroidsZ = finalCentroidsZ;
        end;
        
        if smoothing(1)
            title('Smooth centroid trajectories');
            xlabel('Time point');
            ylabel('Centroid z-position (unit: z-planes)');
            if subtractBaseline
                legend('Raw', 'Baseline-subtracted', 'Baseline-subtracted smooth', 'Baseline-reconstituted smooth (final)', 'Location', 'SouthOutside');
            else
                legend('Raw', 'Smooth (final)', 'Location', 'SouthOutside');
            end;
            h = gca;
            saveas(h, centroidsSmoothedPlotsFileName);
        end;
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'processedMetaData', 'xDriftCorrections', 'yDriftCorrections', 'zDriftCorrections', 'zLowPlaneAdditions', 'zHighPlaneAdditions', ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ');
    end;
elseif smoothing(1)~=0
    if splitXRange
        averagedCentroidsZ1 = 0.5*centroidsZ(:, 1, 1)+0.5*centroidsZ(:, 2, 1);
        averagedCentroidsZ2 = 0.5*centroidsZ(:, 1, 2)+0.5*centroidsZ(:, 2, 2);
        
        if ~isempty(smoothingStart)
            smoothCentroidsZ1 = averagedCentroidsZ1;
            smoothCentroidsZ1(smoothingStart:end) = smooth(averagedCentroidsZ1(smoothingStart:end), smoothing(2), 'rloess');
            smoothCentroidsZ2 = averagedCentroidsZ2;
            smoothCentroidsZ2(smoothingStart:end) = smooth(averagedCentroidsZ2(smoothingStart:end), smoothing(2), 'rloess');
        else
            smoothCentroidsZ1 = smooth(averagedCentroidsZ1, smoothing(2), 'rloess');
            smoothCentroidsZ2 = smooth(averagedCentroidsZ2, smoothing(2), 'rloess');
        end;
        
        centroidsRawAndSmoothedPlotsFileName = ['createMasks.centroidsRawAndSmoothed_SPM' num2str(specimen) outputIdentifier '.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ1, 'Color', [1 0 0], 'LineWidth', 1);
        hold on;
        plot(timepoints, smoothCentroidsZ1, 'Color', [0.5 0 0], 'LineWidth', 1);
        plot(timepoints, averagedCentroidsZ2, 'Color', [0 1 0], 'LineWidth', 1);
        plot(timepoints, smoothCentroidsZ2, 'Color', [0 0.5 0], 'LineWidth', 1);
        title('Raw image centroids and smoothed centroid trajectories');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids (top half)', 'Smooth centroids (top half)', 'Raw centroids (bottom half)', 'Smooth centroids (bottom half)');
        h = gca;
        saveas(h, centroidsRawAndSmoothedPlotsFileName);
        
        if lockLowCentroid(1)
            finalCentroidsZ1 = smoothCentroidsZ1;
            processedCentroidsZ1 = finalCentroidsZ1;
            
            finalCentroidsZ2 = smoothCentroidsZ2;
            finalCentroidsZ2(lockLowCentroid(2):end) = finalCentroidsZ1(lockLowCentroid(2):end) + ...
                finalCentroidsZ2(lockLowCentroid(2)-1) - finalCentroidsZ1(lockLowCentroid(2)-1);
            
            processedCentroidsZ2 = finalCentroidsZ2;
            
            centroidsFixedPlotsFileName = ['createMasks.centroidsFixed_SPM' num2str(specimen) outputIdentifier '.png'];
        
            figure;
            plot(timepoints, smoothCentroidsZ1, 'Color', [0.5 0 0], 'LineWidth', 1);
            hold on;
            plot(timepoints, smoothCentroidsZ2, 'Color', [0 1 0], 'LineWidth', 1);
            plot(timepoints, finalCentroidsZ2, 'Color', [0 0.5 0], 'LineWidth', 1);
            title('Smoothed and fixed centroid trajectories');
            xlabel('Time point');
            ylabel('Z position/shift (unit: z-planes)');
            legend('Smooth centroids (top half, final)', 'Smooth centroids (bottom half)', 'Fixed centroids (bottom half, final)');
            h = gca;
            saveas(h, centroidsFixedPlotsFileName);
        else
            finalCentroidsZ1 = smoothCentroidsZ1;
            processedCentroidsZ1 = finalCentroidsZ1;
            
            finalCentroidsZ2 = smoothCentroidsZ2;
            processedCentroidsZ2 = finalCentroidsZ2;
        end;
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ1', 'processedCentroidsZ2');
    else
        averagedCentroidsZ = 0.5*centroidsZ(:, 1)+0.5*centroidsZ(:, 2);
        
        if ~isempty(smoothingStart)
            smoothCentroidsZ = averagedCentroidsZ;
            smoothCentroidsZ(smoothingStart:end) = smooth(averagedCentroidsZ(smoothingStart:end), smoothing(2), 'rloess');
        else
            smoothCentroidsZ = smooth(averagedCentroidsZ, smoothing(2), 'rloess');
        end;
        
        finalCentroidsZ = smoothCentroidsZ;
        processedCentroidsZ = finalCentroidsZ;
        
        centroidsRawAndSmoothedPlotsFileName = ['createMasks.centroidsRawAndSmoothed_SPM' num2str(specimen) outputIdentifier '.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        hold on;
        plot(timepoints, smoothCentroidsZ, 'k', 'LineWidth', 1);
        title('Raw image centroids and smoothed centroid trajectory');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids', 'Smooth centroids');
        h = gca;
        saveas(h, centroidsRawAndSmoothedPlotsFileName);
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ');
    end;
else
    if splitXRange
        averagedCentroidsZ1 = 0.5*centroidsZ(:, 1, 1)+0.5*centroidsZ(:, 2, 1);
        averagedCentroidsZ2 = 0.5*centroidsZ(:, 1, 2)+0.5*centroidsZ(:, 2, 2);
        
        centroidsRawPlotsFileName = ['createMasks.centroidsRaw_SPM' num2str(specimen) outputIdentifier '.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ1, 'Color', [0.5 0 0], 'LineWidth', 1);
        hold on;
        plot(timepoints, averagedCentroidsZ2, 'Color', [0 0.5 0], 'LineWidth', 1);
        title('Raw image centroids');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids (top half)', 'Raw centroids (bottom half)');
        h = gca;
        saveas(h, centroidsRawPlotsFileName);
        
        if lockLowCentroid(1)
            finalCentroidsZ1 = averagedCentroidsZ1;
            processedCentroidsZ1 = finalCentroidsZ1;
            
            finalCentroidsZ2 = averagedCentroidsZ2;
            finalCentroidsZ2(lockLowCentroid(2):end) = finalCentroidsZ1(lockLowCentroid(2):end) + ...
                finalCentroidsZ2(lockLowCentroid(2)-1) - finalCentroidsZ1(lockLowCentroid(2)-1);
            
            processedCentroidsZ2 = finalCentroidsZ2;
            
            centroidsFixedPlotsFileName = ['createMasks.centroidsFixed_SPM' num2str(specimen) outputIdentifier '.png'];
        
            figure;
            plot(timepoints, averagedCentroidsZ1, 'Color', [0.5 0 0], 'LineWidth', 1);
            hold on;
            plot(timepoints, averagedCentroidsZ2, 'Color', [0 1 0], 'LineWidth', 1);
            plot(timepoints, finalCentroidsZ2, 'Color', [0 0.5 0], 'LineWidth', 1);
            title('Raw and fixed centroid trajectories');
            xlabel('Time point');
            ylabel('Z position/shift (unit: z-planes)');
            legend('Raw centroids (top half, final)', 'Raw centroids (bottom half)', 'Fixed centroids (bottom half, final)');
            h = gca;
            saveas(h, centroidsFixedPlotsFileName);
        else
            finalCentroidsZ1 = averagedCentroidsZ1;
            processedCentroidsZ1 = finalCentroidsZ1;
            
            finalCentroidsZ2 = averagedCentroidsZ2;
            processedCentroidsZ2 = finalCentroidsZ2;
        end;
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ1', 'processedCentroidsZ2');
    else
        averagedCentroidsZ = 0.5*centroidsZ(:, 1)+0.5*centroidsZ(:, 2);
        
        finalCentroidsZ = averagedCentroidsZ;
        processedCentroidsZ = finalCentroidsZ;
        
        centroidsRawPlotsFileName = ['createMasks.centroidsRaw_SPM' num2str(specimen) outputIdentifier '.png'];
        
        figure;
        plot(timepoints, averagedCentroidsZ, 'k', 'LineWidth', 1);
        title('Raw image centroids');
        xlabel('Time point');
        ylabel('Z position/shift (unit: z-planes)');
        legend('Raw centroids');
        h = gca;
        saveas(h, centroidsRawPlotsFileName);
        
        disp(' ');
        disp('Saving global processing logs');
        
        globalLogsFileName = ['createMasks.processingLogs_SPM' num2str(specimen) outputIdentifier '.mat'];
        save(globalLogsFileName, ...
            'xSizes', 'ySizes', 'zSizes', 'minIntensities', 'meanIntensities', 'absoluteThresholds', 'centroidsX', 'centroidsY', 'centroidsZ', ...
            'processedCentroidsZ');
    end;
end;

%% Create MVD masks

disp(' ');
disp('Creating MVD masks');

if transitionType(1) == 0
    transitionVector = (transitionLength:-1:0)./transitionLength;
else
    transitionVector = 1./(1+exp(-((transitionLength/2):-1:(-transitionLength/2)).*10./(transitionType(2)*transitionLength)));
end;
if splitXRange
    transitionZonePositive = repmat(reshape(transitionVector, [1 1 transitionLength+1]), [1 ySize 1]);
    transitionZoneNegative = repmat(reshape(fliplr(transitionVector), [1 1 transitionLength+1]), [1 ySize 1]);
else
    transitionZonePositive = repmat(reshape(transitionVector, [1 1 transitionLength+1]), [xSizes(end, 1) ySize 1]);
    transitionZoneNegative = repmat(reshape(fliplr(transitionVector), [1 1 transitionLength+1]), [xSizes(end, 1) ySize 1]);
end;

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);
disp(' ');

parfor t = 1:numel(timepoints)
    disp(['*** Processing time point ' num2str(timepoints(t))]);
    
    if splitXRange
        averagedCentroidX1 = 0.5*centroidsX(t, 1, 1)+0.5*centroidsX(t, 2, 1);
        averagedCentroidX2 = 0.5*centroidsX(t, 1, 2)+0.5*centroidsX(t, 2, 2);
        if restrictXRange(1)
            currentCenterX1 = averagedCentroidX1+restrictXRange(2)-1;
            currentCenterX2 = averagedCentroidX2+restrictXRange(2)-1;
        else
            xCenter = round(xSizes(numel(timepoints), 1)/2);
            currentCenterX1 = averagedCentroidX1;
            currentCenterX2 = averagedCentroidX2;
        end;
        currentCentroidZ1 = processedCentroidsZ1(t, size(processedCentroidsZ1, 2));
        currentCentroidZ2 = processedCentroidsZ2(t, size(processedCentroidsZ2, 2));
        
        zLocations = ((1:xSizes(numel(timepoints), 1))-currentCenterX1).*((currentCentroidZ2-currentCentroidZ1)/(currentCenterX2-currentCenterX1))+currentCentroidZ1;
    else
        currentCentroidZ = processedCentroidsZ(t, numel(timepoints));
    end;
    
    for camera = [0 1]
        weightsMatrixFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weightsMatrix' outputIdentifier '.klb'];
        weightedStackFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' outputIdentifier '.klb'];
        
        weightedStackXYProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' outputIdentifier '.xyProjection.klb'];
        weightedStackXZProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' outputIdentifier '.xzProjection.klb'];
        weightedStackYZProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' outputIdentifier '.yzProjection.klb'];
        
        if exist(weightsMatrixFilename, 'file') ~= 2 || (applyWeights && exist(weightedStackFilename, 'file') ~= 2) || (applyWeights && saveProjections && exist(weightedStackYZProjectionFilename, 'file') ~= 2) || forceWeights
            disp(['    Creating and/or applying weights matrix for camera ' num2str(camera)]);
            
            stackFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
                filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') inputIdentifier '.klb'];
            stack = readImage(stackFilename);
            
            if exist(weightsMatrixFilename, 'file') ~= 2 || forceWeights
                weights = zeros(xSizes(t, camera+1), ySizes(t, camera+1), zSizes(t, camera+1), 'single');
                
                if camera == primaryCamera
                    if splitXRange
                        zLocationsInteger = round(zLocations);
                        zLocationsStepIndices = unique([1 find(zLocationsInteger(1:(end-1)) ~= zLocationsInteger(2:end))+1 xSizes(end, 1)]);
                        
                        weightsXZAssignments = zeros(xSizes(t, camera+1), zSizes(t, camera+1), 'uint8');
                        weightsXZSegment = zeros(xSizes(t, camera+1), zSizes(t, camera+1), 'single');
                        
                        for x = 1:xSizes(end, 1)
                            currentStartZ = round(zLocations(x))+transitionDelta;
                            currentStopZ = round(zLocations(x))+transitionDelta+transitionLength;
                            if currentStartZ<2
                                currentStartZ = 2;
                                currentStopZ = 2+transitionLength;
                            elseif currentStopZ>(zSizes(t, camera+1)-1)
                                currentStartZ = zSizes(t, camera+1)-1-transitionLength;
                                currentStopZ = zSizes(t, camera+1)-1;
                            end;
                            
                            weightsXZSegment(x, 1:(currentStartZ-1)) = 1;
                            weightsXZAssignments(x, 1:(currentStartZ-1)) = 1;
                            weightsXZAssignments(x, (currentStopZ+1):end) = 1;
                            if ~isempty(find(zLocationsStepIndices == x, 1))
                                weightsXZSegment(x, currentStartZ:currentStopZ) = transitionZonePositive(1, 1, :);
                                weightsXZAssignments(x, currentStartZ:currentStopZ) = 1;
                            else
                                weightsXZAssignments(x, currentStartZ:currentStopZ) = 2;
                            end;
                        end;
                        
                        for z = (min(zLocationsInteger)+transitionDelta):(max(zLocationsInteger)+transitionDelta+transitionLength)
                            sliverAssignments = weightsXZAssignments(:, z);
                            sliverSegment = weightsXZSegment(:, z);
                            if ~isempty(find(sliverAssignments == 2, 1))
                                firstSlot = find(sliverAssignments == 2, 1, 'first');
                                lastSlot = find(sliverAssignments == 2, 1, 'last');
                                if sliverSegment(firstSlot-1) == sliverSegment(lastSlot+1)
                                    weightsXZSegment((firstSlot-1):(lastSlot+1), z) = sliverSegment(firstSlot-1);
                                else
                                    weightsXZSegment((firstSlot-1):(lastSlot+1), z) = ...
                                        (sliverSegment(firstSlot-1):((sliverSegment(lastSlot+1)-sliverSegment(firstSlot-1))/(lastSlot-firstSlot+2)):sliverSegment(lastSlot+1))';
                                end;
                            end;
                        end;
                        
                        for y = 1:ySizes(t, camera+1)
                            weights(:, y, :) = weightsXZSegment;
                        end;
                    else
                        currentStartZ = round(currentCentroidZ)+transitionDelta;
                        currentStopZ = round(currentCentroidZ)+transitionDelta+transitionLength;
                        
                        if currentStartZ<1
                            currentStartZ = 1;
                            currentStopZ = 1+transitionLength;
                        elseif currentStopZ>zSizes(t, camera+1)
                            currentStartZ = zSizes(t, camera+1)-transitionLength;
                            currentStopZ = zSizes(t, camera+1);
                        end;
                        
                        weights(:, :, 1:(currentStartZ-1)) = 1;
                        weights(:, :, currentStartZ:currentStopZ) = transitionZonePositive;
                    end;
                else
                    if splitXRange
                        zLocationsInteger = round(zLocations);
                        zLocationsStepIndices = unique([1 find(zLocationsInteger(1:(end-1)) ~= zLocationsInteger(2:end))+1 xSizes(end, 1)]);
                        
                        weightsXZAssignments = zeros(xSizes(t, camera+1), zSizes(t, camera+1), 'uint8');
                        weightsXZSegment = zeros(xSizes(t, camera+1), zSizes(t, camera+1), 'single');
                        
                        for x = 1:xSizes(end, 1)
                            currentStartZ = round(zLocations(x))-transitionDelta;
                            currentStopZ = round(zLocations(x))-transitionDelta-transitionLength;
                            if currentStopZ<2
                                currentStartZ = 2+transitionLength;
                                currentStopZ = 2;
                            elseif currentStartZ>(zSizes(t, camera+1)-1)
                                currentStartZ = zSizes(t, camera+1)-1;
                                currentStopZ = zSizes(t, camera+1)-1-transitionLength;
                            end;
                            
                            weightsXZSegment(x, (currentStartZ+1):end) = 1;
                            weightsXZAssignments(x, 1:(currentStopZ-1)) = 1;
                            weightsXZAssignments(x, (currentStartZ+1):end) = 1;
                            if ~isempty(find(zLocationsStepIndices == x, 1))
                                weightsXZSegment(x, currentStopZ:currentStartZ) = transitionZoneNegative(1, 1, :);
                                weightsXZAssignments(x, currentStopZ:currentStartZ) = 1;
                            else
                                weightsXZAssignments(x, currentStopZ:currentStartZ) = 2;
                            end;
                        end;
                        
                        for z = (min(zLocationsInteger)-transitionDelta-transitionLength):(max(zLocationsInteger)-transitionDelta)
                            sliverAssignments = weightsXZAssignments(:, z);
                            sliverSegment = weightsXZSegment(:, z);
                            if ~isempty(find(sliverAssignments == 2, 1))
                                firstSlot = find(sliverAssignments == 2, 1, 'first');
                                lastSlot = find(sliverAssignments == 2, 1, 'last');
                                if sliverSegment(firstSlot-1) == sliverSegment(lastSlot+1)
                                    weightsXZSegment((firstSlot-1):(lastSlot+1), z) = sliverSegment(firstSlot-1);
                                else
                                    weightsXZSegment((firstSlot-1):(lastSlot+1), z) = ...
                                        (sliverSegment(firstSlot-1):((sliverSegment(lastSlot+1)-sliverSegment(firstSlot-1))/(lastSlot-firstSlot+2)):sliverSegment(lastSlot+1))';
                                end;
                            end;
                        end;
                        
                        for y = 1:ySizes(t, camera+1)
                            weights(:, y, :) = weightsXZSegment;
                        end;
                    else
                        currentStartZ = round(currentCentroidZ)-transitionDelta;
                        currentStopZ = round(currentCentroidZ)-transitionDelta-transitionLength;
                        
                        if currentStopZ<1
                            currentStartZ = 1+transitionLength;
                            currentStopZ = 1;
                        elseif currentStartZ>zSizes(t, camera+1)
                            currentStartZ = zSizes(t, camera+1);
                            currentStopZ = zSizes(t, camera+1)-transitionLength;
                        end;
                        
                        weights(:, :, (currentStopZ+1):end) = 1;
                        weights(:, :, currentStopZ:currentStartZ) = transitionZoneNegative;
                    end;
                end;
                
                writeImage(weights, weightsMatrixFilename);
            else
                disp(['    Loading weights matrix for camera ' num2str(camera)]);
                
                weights = readImage(weightsMatrixFilename);
            end;
            
            if applyWeights
                if transitionToZero(1)
                    weightedStack = uint16(single(stack).*weights);
                else
                    if transitionToZero(2) == 0
                        currentBackground = minIntensities(t, camera+1);
                    else
                        currentBackground = transitionToZero(2);
                    end;
                    weightedStack = uint16(single(stack-currentBackground).*weights + currentBackground);
                end;
                % clear stack weights;
                writeImage(weightedStack, weightedStackFilename);
                
                if saveProjections
                    xyProjection = max(weightedStack, [], 3);
                    xzProjection = squeeze(max(weightedStack, [], 2));
                    yzProjection = squeeze(max(weightedStack, [], 1));
                    % clear weightedStack;
                    
                    writeImage(xyProjection, weightedStackXYProjectionFilename);
                    writeImage(xzProjection, weightedStackXZProjectionFilename);
                    writeImage(yzProjection, weightedStackYZProjectionFilename);
                end;
            end;
        else
            disp(['    Existing results found for camera ' num2str(camera)]);
        end;
    end;
end;

disp(' ');
if matlabpool('size') > 0
    matlabpool('close');
end;
disp(' ');