function timeFuse(parameterDatabase, t, memoryEstimate)

% -----------------------------------------------------------------------------------------------
% | Time-lapse image data fusion for simultaneous multi-view light sheet microscopy             |
% |                                                                                             |
% | Bi-directional illumination dual-camera image stack registration and fusion,                |
% | using adaptive slicing, correlation by rigid transformation and adaptive fusion             |
% | by linear short-distance blending, arithmetic averaging or wavelet decomposition            |
% |                                                                                             |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2015                          |
% | Email: kellerp@janelia.hhmi.org                                                             |
% |                                                                                             |
% | Utilizes optimization modules and functions by Fernando Amat, HHMI/Janelia Research Campus: |
% | imgaussianAnisotropy.m (included below)                                                     |
% | fminuncFA.m            (included below)                                                     |
% | transformCamera.m      (included below)                                                     |
% | transformChannel.m     (included below)                                                     |
% | readKLBstack.mexw64                                                                         |
% | writeKLBstack.mexw64                                                                        |
% |                                                                                             |
% | Utilizes optimization modules from the MathWorks File Exchange:                             |
% | fInterpolate.mexw64 (ba_interp2.mexw64)                                                     |
% -----------------------------------------------------------------------------------------------

% lookUpTable format (4-view fusion): (will be reduced to rows 2-14)
% row 1: time point
% row 2: camera 0 bestOffset
% row 3: camera 0 bestRotation
% row 4: camera 0 correctionFactor
% row 5: camera 0 correctionFlag
% row 6: camera 1 bestOffset
% row 7: camera 1 bestRotation
% row 8: camera 1 correctionFactor
% row 9: camera 1 correctionFlag
% row 10: final bestXOffset
% row 11: final bestYOffset
% row 12: final bestRotation
% row 13: final correctionFactor
% row 14: final correctionFlag

% lookUpTable format (2-view channel fusion): (will be reduced to rows 2-5)
% row 1: time point
% row 2: camera X bestOffset
% row 3: camera X bestRotation
% row 4: camera X correctionFactor
% row 5: camera X correctionFlag

% lookUpTable format (2-view camera fusion): (will be reduced to rows 2-6)
% row 1: time point
% row 2: final bestXOffset
% row 3: final bestYOffset
% row 4: final bestRotation
% row 5: final correctionFactor
% row 6: final correctionFlag

load(parameterDatabase);
timepoint = timepoints(t);

load(lookUpTable);

version = 1.17;

configuration = cell(33, 1);

configuration{1}  = version;      configuration{2}  = timepoint;    configuration{3}  = references;   configuration{4}  = globalMask;
configuration{5}  = inputString;  configuration{6}  = sourceString; configuration{7}  = outputString;
configuration{8}  = inputID;      configuration{9}  = outputID;     configuration{10} = lookUpTable;  configuration{11} = dataType;
configuration{12} = specimen;     configuration{13} = cameras;      configuration{14} = sChannels;    configuration{15} = tChannels;
configuration{16} = reducedIO;    configuration{17} = inputType;    configuration{18} = outputType;
configuration{19} = splitting;    configuration{20} = intSizes;     configuration{21} = correction;
configuration{22} = percentile;   configuration{23} = subSampling;  configuration{24} = fusionType;   configuration{25} = blending;
configuration{26} = enforceFlag;  configuration{27} = cropping;     configuration{28} = scaling;
configuration{29} = leftFlags;    configuration{30} = flipHFlag;    configuration{31} = flipVFlag;    configuration{32} = frontFlag;
configuration{33} = [jobMemory(1) memoryEstimate];

%% main loop

if length(cameras) == 2 && length(tChannels) == 2
    processingMode = 0; % 4-view fusion
elseif length(cameras) == 1 && length(tChannels) == 2
    processingMode = 1; % 2-view channel fusion
elseif length(cameras) == 2 && length(tChannels) == 1
    processingMode = 2; % 2-view camera fusion
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion)');
end;

if verbose == 1 || verbose == 2
    disp(' ');
end;

jobCompleted = 0;

disp(['processing time point ' num2str(timepoint, '%.6d')]);

distanceArray = abs(references - timepoint);
currentSource = references(find(distanceArray == min(distanceArray), 1));

inputFolder   = [inputString '/SPM' num2str(specimen, '%.2d') '/TM' num2str(timepoint, '%.6d')];
inputHeader   = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d')];
sourceFolder  = [sourceString '.TM' num2str(currentSource, '%.6d') '_multiFused' inputID];
globalFolder  = [sourceString '.TM' num2str(globalMask(2), '%.6d') '_multiFused' inputID];
outputFolder  = [outputString '.TM' num2str(timepoint, '%.6d') '_timeFused' outputID];
outputHeader  = [outputFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d')];

currentTransformations = lookUpTable(find(lookUpTable(:, 1) == timepoint, 1), 2:end);

if verbose == 1 || verbose == 2
    disp(' ');
end;

if ~isdir(outputFolder)
    mkdir(outputFolder);
end;

switch processingMode
    case 0
        outputString = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d')];
        primaryDataArray = cell(2, 2);
    case 1
        outputString = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d')];
        primaryDataArray = cell(1, 2);
    case 2
        outputString = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d')];
        primaryDataArray = cell(2, 1);
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

save([outputString '.configuration.mat'], 'configuration');

if processingMode ~= 2
    
    %% alignChannels.m
    
    tic;
    
    for camera = cameras
        currentCamera = find(cameras == camera, 1);
        
        sX = cropping{find(cameras == camera, 1)}(2) + 1;
        eX = cropping{find(cameras == camera, 1)}(2) + cropping{find(cameras == camera, 1)}(4);
        sY = cropping{find(cameras == camera, 1)}(1) + 1;
        eY = cropping{find(cameras == camera, 1)}(1) + cropping{find(cameras == camera, 1)}(3);
        sZ = cropping{find(cameras == camera, 1)}(5) + 1;
        eZ = cropping{find(cameras == camera, 1)}(6) + 1;
        
        xSize = eX - sX + 1;
        ySize = eY - sY + 1;
        zSize = eZ - sZ + 1;
        
        for c = tChannels
            if verbose == 2
                disp('--');
                disp(['processing camera ' num2str(camera) ' channel ' num2str(c)]);
                disp('--');
            end;
            
            currentChannel = find(tChannels == c, 1);
            
            header = [inputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(c, '%.2d')];            
            stackName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(c, '%.2d') '.stack' outputExtension];
            
            if verbose == 2
                disp('reading image stack');
            end;
            
            if exist(stackName, 'file') ~= 2 || enforceFlag(1)
                primaryDataArray{currentCamera, currentChannel} = readImage([inputFolder '/' header inputExtension]);
                if ~strcmp(class(primaryDataArray{currentCamera, currentChannel}), 'uint16')
                    primaryDataArray{currentCamera, currentChannel} = uint16(primaryDataArray{currentCamera, currentChannel});
                end;
                
                if sum(cropping{find(cameras == camera, 1)}) > 0
                    primaryDataArray{currentCamera, currentChannel} = primaryDataArray{currentCamera, currentChannel}(sX:eX, sY:eY, sZ:eZ);
                else
                    xSize = size(primaryDataArray{currentCamera, currentChannel}, 1);
                    ySize = size(primaryDataArray{currentCamera, currentChannel}, 2);
                    zSize = size(primaryDataArray{currentCamera, currentChannel}, 3);
                end;
                
                if ~reducedIO
                    writeImage(primaryDataArray{currentCamera, currentChannel}, stackName);
                end;
            else
                primaryDataArray{currentCamera, currentChannel} = readImage(stackName);
                
                if sum(cropping{find(cameras == camera, 1)}) == 0
                    xSize = size(primaryDataArray{currentCamera, currentChannel}, 1);
                    ySize = size(primaryDataArray{currentCamera, currentChannel}, 2);
                    zSize = size(primaryDataArray{currentCamera, currentChannel}, 3);
                end;
            end;
            
            if verbose == 2
                disp(' ');
            end;
        end;
        
        if verbose == 2
            disp('--');
            disp(['transforming camera ' num2str(camera) ' data']);
            disp('--');
        end;
        
        transformationIndex = 1 + 4 * (find(cameras == camera, 1) - 1);
        bestOffset   = currentTransformations(transformationIndex);
        bestRotation = currentTransformations(transformationIndex + 1);
        
        transformedStackName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.transformedStack' outputExtension];
        
        if exist(transformedStackName, 'file') ~= 2 || enforceFlag(2)
            if verbose == 2
                disp('transforming image stack');
            end;
            
            RR = [ cosd(bestRotation) sind(bestRotation);...
                -sind(bestRotation) cosd(bestRotation)];
            [XI ZI] = meshgrid(1:xSize, 1:zSize);
            cc = [(xSize + 1) / 2.0 (zSize + 1) / 2.0];
            XI = XI' - cc(1);
            ZI = ZI' - cc(2);
            XZaux = (RR * [XI(:) ZI(:) .* scaling]')';
            XI = reshape(XZaux(:, 1) + cc(1), size(XI));
            ZI = reshape(XZaux(:, 2) ./ scaling + cc(2), size(ZI));
            
            ZI = ZI - bestOffset;
            
            maskPos = find(XI < 1 | ZI < 1 | XI > xSize | ZI > zSize);
            maskTR = true(xSize, zSize);
            maskTR(maskPos) = false;
            
            permutedStack = permute(primaryDataArray{currentCamera, 2}, [1 3 2]);
            
            primaryDataArray{currentCamera, 2} = zeros(xSize, zSize, ySize, 'uint16');
            for i = 1:splitting
                ySlabStart = round((i - 1) * ySize / splitting + 1);
                ySlabStop  = round(i * ySize / splitting);
                primaryDataArray{currentCamera, 2}(:, :, ySlabStart:ySlabStop) = ...
                    uint16(bsxfun(@times, maskTR, fInterpolate(double(permutedStack(:, :, ySlabStart:ySlabStop)), ZI, XI, 'cubic')));
            end;
            
            primaryDataArray{currentCamera, 2} = permute(primaryDataArray{currentCamera, 2}, [1 3 2]);
            
            clear permutedStack;
            
            if ~reducedIO
                writeImage(primaryDataArray{currentCamera, 2}, transformedStackName);
            end;
        else
            if verbose == 2
                disp('reading transformed stack');
            end;
            
            primaryDataArray{currentCamera, 2} = readImage(transformedStackName);
        end;
        
        if verbose == 2
            disp(' ');
        end;
    end;
    
    time = toc;
    if verbose == 1 || verbose == 2
        disp(['*** channel alignment completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
    end;
    
    if verbose == 2
        disp(' ');
    end;
    
    %% fuseChannels.m
    
    tic;
    
    for camera = cameras
        if verbose == 2
            disp('--');
            disp(['processing camera ' num2str(camera) ' data']);
            disp('--');
        end;
        
        currentCamera = find(cameras == camera, 1);
        
        leftFlag = leftFlags(find(cameras == camera, 1));
        
        xSize = size(primaryDataArray{currentCamera, 1}, 1);
        ySize = size(primaryDataArray{currentCamera, 1}, 2);
        zSize = size(primaryDataArray{currentCamera, 1}, 3);
        
        if verbose == 2
            disp('reading average mask');
        end;
        
        if globalMask(1) == 0
            averageMaskName = [sourceFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(currentSource, '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '_CHN' num2str(sChannels(2), '%.2d') '.fusionMask' inputExtension];
        else
            averageMaskName = [globalFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(globalMask(2), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '_CHN' num2str(sChannels(2), '%.2d') '.fusionMask' inputExtension];
        end;
        averageMask = readImage(averageMaskName);
        
        if correction(1) ~= 0
            if dataType == 1
                minIntensityName = [inputFolder '/' inputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.minIntensity.mat'];
                load(minIntensityName, 'minIntensity');
                backgroundIntensity1 = minIntensity(end);
                
                minIntensityName = [inputFolder '/' inputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.minIntensity.mat'];
                load(minIntensityName, 'minIntensity');
                backgroundIntensity2 = minIntensity(end);
            else
                backgroundArray1 = primaryDataArray{currentCamera, 1}(1:subSampling:end);
                backgroundArray2 = primaryDataArray{currentCamera, 2}(1:subSampling:end);
                background = prctile(cat(2, backgroundArray1(backgroundArray1 > 0), backgroundArray2(backgroundArray2 > 0)), percentile);
                clear backgroundArray1 backgroundArray2;
            end;
            
            if correction(1) == 1
                if verbose == 2
                    disp('analyzing intensity difference');
                end;
                
                minMask = min(averageMask(averageMask > 0));
                maxMask = max(averageMask(:));
                
                lowerBound = double(minMask) - double(intSizes(1));
                if lowerBound < 1
                    lowerSlabSize = intSizes(1) + lowerBound - 1;
                else
                    lowerSlabSize = intSizes(1);
                end;
                
                upperBound = double(maxMask) + double(intSizes(1)) - double(ySize);
                if upperBound > 0
                    upperSlabSize = intSizes(1) - upperBound;
                else
                    upperSlabSize = intSizes(1);
                end;
                
                finalSlabSize = lowerSlabSize + upperSlabSize + 1;
                
                dataSlice1 = zeros(xSize, zSize, finalSlabSize, 'uint16');
                dataSlice2 = zeros(xSize, zSize, finalSlabSize, 'uint16');
                for z = 1:zSize
                    for x = 1:xSize
                        if averageMask(x, z) > 0
                            currentIndexArray = (averageMask(x, z) - lowerSlabSize):(averageMask(x, z) + upperSlabSize);
                            dataSlice1(x, z, :) = primaryDataArray{currentCamera, 1}(x, currentIndexArray, z);
                            dataSlice2(x, z, :) = primaryDataArray{currentCamera, 2}(x, currentIndexArray, z);
                        end;
                    end;
                end;
                
                if dataType == 1
                    dataSlice1 = dataSlice1 - backgroundIntensity1;
                    dataSlice2 = dataSlice2 - backgroundIntensity2;
                else
                    dataSlice1 = dataSlice1 - background;
                    dataSlice2 = dataSlice2 - background;
                end;
                
                dataSlice1Name = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.fusionDataSlice' outputExtension];
                dataSlice2Name = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.fusionDataSlice' outputExtension];
                
                if ~reducedIO || reducedIO == 2
                    writeImage(dataSlice1, dataSlice1Name);
                    writeImage(dataSlice2, dataSlice2Name);
                end;
                
                intensityList1 = double(dataSlice1(:));
                intensitySum1 = sum(intensityList1(intensityList1 > 0));
                intensityList2 = double(dataSlice2(:));
                intensitySum2 = sum(intensityList2(intensityList2 > 0));
                
                if intensitySum1 > intensitySum2
                    correctionFactor = intensitySum1 / intensitySum2;
                    correctionFlag = 2;
                else
                    correctionFactor = intensitySum2 / intensitySum1;
                    correctionFlag = 1;
                end;
                
                if dataType == 1
                    intensityCorrection = [double(backgroundIntensity1) double(backgroundIntensity2) double(intensitySum1) double(intensitySum2) double(correctionFactor) double(correctionFlag)];
                else
                    intensityCorrection = [double(background) double(intensitySum1) double(intensitySum2) double(correctionFactor) double(correctionFlag)];
                end;
                intensityCorrectionName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.intensityCorrection.mat'];
                save(intensityCorrectionName, 'intensityCorrection');
            elseif correction(1) == 2
                if correction(3) > 1
                    correctionFactor = correction(3);
                    correctionFlag = 2;
                else
                    correctionFactor = 1 / correction(3);
                    correctionFlag = 1;
                end;
            elseif correction(1) == 3
                transformationIndex = 3 + 4 * (find(cameras == camera, 1) - 1);
                correctionFactor    = currentTransformations(transformationIndex);
                correctionFlag      = currentTransformations(transformationIndex + 1);
            end;
            
            correctedStack1Name = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.correctedStack' outputExtension];
            correctedStack2Name = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.correctedStack' outputExtension];
            
            if exist(correctedStack1Name, 'file') ~= 2 || exist(correctedStack2Name, 'file') ~= 2 || enforceFlag(3)
                if verbose == 2
                    disp('correcting interpolated image stacks');
                end;
                
                if dataType == 1
                    primaryDataArray{currentCamera, 1} = primaryDataArray{currentCamera, 1} - backgroundIntensity1;
                    primaryDataArray{currentCamera, 1}(primaryDataArray{currentCamera, 1} < 0) = 0;
                    primaryDataArray{currentCamera, 2} = primaryDataArray{currentCamera, 2} - backgroundIntensity2;
                    primaryDataArray{currentCamera, 2}(primaryDataArray{currentCamera, 2} < 0) = 0;
                else
                    primaryDataArray{currentCamera, 1} = primaryDataArray{currentCamera, 1} - background;
                    primaryDataArray{currentCamera, 1}(primaryDataArray{currentCamera, 1} < 0) = 0;
                    primaryDataArray{currentCamera, 2} = primaryDataArray{currentCamera, 2} - background;
                    primaryDataArray{currentCamera, 2}(primaryDataArray{currentCamera, 2} < 0) = 0;
                end;
                
                if correctionFlag == 1
                    primaryDataArray{currentCamera, 1} = primaryDataArray{currentCamera, 1} .* correctionFactor;
                else
                    primaryDataArray{currentCamera, 2} = primaryDataArray{currentCamera, 2} .* correctionFactor;
                end;
                
                if ~reducedIO
                    writeImage(primaryDataArray{currentCamera, 1}, correctedStack1Name);
                    writeImage(primaryDataArray{currentCamera, 2}, correctedStack2Name);
                end;
            else
                if verbose == 2
                    disp('reading intensity-corrected stacks');
                end;
                
                primaryDataArray{currentCamera, 1} = readImage(correctedStack1Name);
                primaryDataArray{currentCamera, 2} = readImage(correctedStack2Name);
            end;
        else
            if verbose == 2
                disp('skipping intensity correction');
            end;
        end;
        
        if fusionType == 0
            stitchedStackName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.stitchedStack' outputExtension];
            
            if exist(stitchedStackName, 'file') ~= 2 || enforceFlag(4)
                if verbose == 2
                    disp('stitching corrected stacks');
                end;
                
                channelMask = bsxfun(@le, reshape(uint16(1:ySize), [1, ySize, 1]), reshape(averageMask, [xSize, 1, zSize]));
                if leftFlag == 1
                    fusedStack = bsxfun(@times, uint16(reshape(averageMask, [xSize 1 zSize]) > 0), primaryDataArray{currentCamera, 2});
                    fusedStack(channelMask) = primaryDataArray{currentCamera, 1}(channelMask);
                else
                    fusedStack = bsxfun(@times, uint16(reshape(averageMask, [xSize 1 zSize]) > 0), primaryDataArray{currentCamera, 1});
                    fusedStack(channelMask) = primaryDataArray{currentCamera, 2}(channelMask);
                end;
                clear channelMask;
                
                if ~reducedIO
                    writeImage(fusedStack, stitchedStackName);
                end;
            else
                if verbose == 2
                    disp('reading stitched stack');
                end;
                
                fusedStack = readImage(stitchedStackName);
            end;
        end;
        
        fusedStackName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.fusedStack' outputExtension];
        
        if exist(fusedStackName, 'file') ~= 2 || enforceFlag(5)
            if fusionType == 0
                if verbose == 2
                    disp('blending corrected stacks');
                end;
                
                weighting1Left = 1:((0.5 - 1) / (blending(1) - 1)):0.5; % dominant fraction
                weighting2Left = 0:((0.5 - 0) / (blending(1) - 1)):0.5; % fading fraction
                
                weighting1Right = 0.5:((1 - 0.5) / (blending(1) - 1)):1; % dominant fraction
                weighting2Right = 0.5:((0 - 0.5) / (blending(1) - 1)):0; % fading fraction
                
                if leftFlag == 1
                    for z = 1:zSize
                        for x = 1:xSize
                            if averageMask(x, z) > 0
                                array1Left = double(primaryDataArray{currentCamera, 1}(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z));
                                array2Left = double(primaryDataArray{currentCamera, 2}(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z));
                                array1Right = double(primaryDataArray{currentCamera, 1}(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z));
                                array2Right = double(primaryDataArray{currentCamera, 2}(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z));
                                fusedStack(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z) = uint16(array1Left .* weighting1Left + array2Left .* weighting2Left);
                                fusedStack(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z) = uint16(array1Right .* weighting2Right + array2Right .* weighting1Right);
                            end;
                        end;
                    end;
                else
                    for z = 1:zSize
                        for x = 1:xSize
                            if averageMask(x, z) > 0
                                array1Left = double(primaryDataArray{currentCamera, 1}(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z));
                                array2Left = double(primaryDataArray{currentCamera, 2}(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z));
                                array1Right = double(primaryDataArray{currentCamera, 1}(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z));
                                array2Right = double(primaryDataArray{currentCamera, 2}(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z));
                                fusedStack(x, (averageMask(x, z) - blending(1) + 1):averageMask(x, z), z) = uint16(array1Left .* weighting2Left + array2Left .* weighting1Left);
                                fusedStack(x, (averageMask(x, z) + 1):(averageMask(x, z) + blending(1)), z) = uint16(array1Right .* weighting1Right + array2Right .* weighting2Right);
                            end;
                        end;
                    end;
                end;
            elseif fusionType == 1
                if verbose == 2
                    disp('performing wavelet fusion of corrected stacks');
                end;
                
                fusedStack = zeros(xSize, ySize, zSize, 'uint16');
                for z = 1:zSize
                    fusedStack(:, :, z) = uint16(wfusimg(primaryDataArray{currentCamera, 1}(:, :, z), primaryDataArray{currentCamera, 2}(:, :, z), 'db4', 5, 'mean', 'max'));
                end;
            elseif fusionType == 2
                if verbose == 2
                    disp('performing averaging of corrected stacks');
                end;
                
                fusedStack = zeros(xSize, ySize, zSize, 'uint16');
                for z = 1:zSize
                    fusedStack(:, :, z) = (primaryDataArray{currentCamera, 1}(:, :, z) + primaryDataArray{currentCamera, 2}(:, :, z)) ./ 2;
                end;
            end;
            
            if (processingMode == 0 && ~reducedIO) || processingMode == 1
                writeImage(fusedStack, fusedStackName);
            end;
            
            if processingMode == 1
                currentString = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d')];
                
                fusedStackProjXYName = [currentString '.fusedStack_xyProjection' outputExtension];
                fusedStackProjXZName = [currentString '.fusedStack_xzProjection' outputExtension];
                fusedStackProjYZName = [currentString '.fusedStack_yzProjection' outputExtension];
                
                writeImage(max(fusedStack, [], 3), fusedStackProjXYName);
                writeImage(squeeze(max(fusedStack, [], 2)), fusedStackProjXZName);
                writeImage(squeeze(max(fusedStack, [], 1)), fusedStackProjYZName);
            end;
        else
            if verbose == 2
                disp('reading fused stack');
            end;
            
            fusedStack = readImage(fusedStackName);
        end;
        
        primaryDataArray{currentCamera, 1} = fusedStack;
        primaryDataArray{currentCamera, 2} = [];
        
        clear fusedStack;
        
        if verbose == 2
            disp(' ');
        end;
    end;
    
    time = toc;
    if verbose == 1
        disp(['***    channel fusion completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
    elseif verbose == 2
        disp(['*** channel fusion completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
    end;
    
    if verbose == 2
        disp(' ');
    end;
end;

if processingMode ~= 1
    
    %% alignCameras.m
    
    tic;
    
    if processingMode == 2
        for camera = cameras
            currentCamera = find(cameras == camera, 1);
            
            sX = cropping{find(cameras == camera, 1)}(2) + 1;
            eX = cropping{find(cameras == camera, 1)}(2) + cropping{find(cameras == camera, 1)}(4);
            sY = cropping{find(cameras == camera, 1)}(1) + 1;
            eY = cropping{find(cameras == camera, 1)}(1) + cropping{find(cameras == camera, 1)}(3);
            sZ = cropping{find(cameras == camera, 1)}(5) + 1;
            eZ = cropping{find(cameras == camera, 1)}(6) + 1;
            
            xSize = eX - sX + 1;
            ySize = eY - sY + 1;
            zSize = eZ - sZ + 1;
            
            if verbose == 2
                disp('--');
                disp(['processing camera ' num2str(camera) ' channel ' num2str(tChannels(1))]);
                disp('--');
            end;
            
            header = [inputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d')];
            stackName = [outputHeader '_CM' num2str(camera, '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.stack' outputExtension];
            
            if verbose == 2
                disp('reading image stack');
            end;
            
            if exist(stackName, 'file') ~= 2 || enforceFlag(1)
                primaryDataArray{currentCamera, 1} = readImage([inputFolder '/' header inputExtension]);
                if ~strcmp(class(primaryDataArray{currentCamera, 1}), 'uint16')
                    primaryDataArray{currentCamera, 1} = uint16(primaryDataArray{currentCamera, 1});
                end;
                
                if sum(cropping{find(cameras == camera, 1)}) > 0
                    primaryDataArray{currentCamera, 1} = primaryDataArray{currentCamera, 1}(sX:eX, sY:eY, sZ:eZ);
                else
                    xSize = size(primaryDataArray{currentCamera, 1}, 1);
                    ySize = size(primaryDataArray{currentCamera, 1}, 2);
                    zSize = size(primaryDataArray{currentCamera, 1}, 3);
                end;
                
                if ~reducedIO
                    writeImage(primaryDataArray{currentCamera, 1}, stackName);
                end;
            else
                primaryDataArray{currentCamera, 1} = readImage(stackName);
                
                if sum(cropping{find(cameras == camera, 1)}) == 0
                    xSize = size(primaryDataArray{currentCamera, 1}, 1);
                    ySize = size(primaryDataArray{currentCamera, 1}, 2);
                    zSize = size(primaryDataArray{currentCamera, 1}, 3);
                end;
            end;
            
            if verbose == 2
                disp(' ');
            end;
        end;
    end;
    
    if flipHFlag
        for z = 1:zSize
            primaryDataArray{2, 1}(:, :, z) = fliplr(primaryDataArray{2, 1}(:, :, z));
        end;
    end;
    
    if flipVFlag
        for z = 1:zSize
            primaryDataArray{2, 1}(:, :, z) = flipud(primaryDataArray{2, 1}(:, :, z));
        end;
    end;
    
    if processingMode == 0
        bestXOffset  = currentTransformations(9);
        bestYOffset  = currentTransformations(10);
        bestRotation = currentTransformations(11);
    else
        bestXOffset  = currentTransformations(1);
        bestYOffset  = currentTransformations(2);
        bestRotation = currentTransformations(3);
    end;
    
    if verbose == 2
        disp('--');
        disp(['transforming camera ' num2str(cameras(2)) ' data']);
        disp('--');
    end;
    
    if processingMode == 0
        transformedStackName = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.transformedStack' outputExtension];
    else
        transformedStackName = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.transformedStack' outputExtension];
    end;
    
    if exist(transformedStackName, 'file') ~= 2 || enforceFlag(6)
        if verbose == 2
            disp('transforming image stack');
        end;
        
        stackSlice = primaryDataArray{2, 1}(:, :, 1);
        
        [XI YI] = meshgrid(1:size(stackSlice, 1), 1:size(stackSlice, 2));
        
        RR = [ cosd(bestRotation) sind(bestRotation);...
            -sind(bestRotation) cosd(bestRotation)];
        cc = [(size(stackSlice, 1) + 1) / 2.0 (size(stackSlice, 2) + 1) / 2.0];
        XI = XI' - cc(1);
        YI = YI' - cc(2);
        XZaux = (RR * [XI(:) YI(:)]')';
        XI = reshape(XZaux(:, 1) + cc(1), size(XI));
        YI = reshape(XZaux(:, 2) + cc(2), size(YI));
        
        XI = XI - bestXOffset;
        YI = YI - bestYOffset;
        
        maskPos = find(XI < 1 | YI < 1 | XI > size(stackSlice, 1) | YI > size(stackSlice, 2));
        maskTR = true(size(stackSlice, 1), size(stackSlice, 2));
        maskTR(maskPos) = false;
        
        transformedStack = zeros(xSize, ySize, zSize, 'uint16');
        for i = 1:splitting
            zSlabStart = round((i - 1) * zSize / splitting + 1);
            zSlabStop  = round(i * zSize / splitting);
            transformedStack(:, :, zSlabStart:zSlabStop) = ...
                uint16(bsxfun(@times, maskTR, fInterpolate(double(primaryDataArray{2, 1}(:, :, zSlabStart:zSlabStop)), YI, XI, 'cubic')));
        end;
        
        if ~reducedIO
            writeImage(transformedStack, transformedStackName);
        end;
    else
        if verbose == 2
            disp('reading transformed image stack');
        end;
        
        transformedStack = readImage(transformedStackName);
    end;
    
    primaryDataArray{2, 1} = transformedStack;
    clear transformedStack;
    
    if verbose == 2
        disp(' ');
    end;
    
    time = toc;
    if verbose == 1
        disp(['***  camera alignment completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
    elseif verbose == 2
        disp(['*** camera alignment completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
    end;
    
    if verbose == 2
        disp(' ');
    end;
    
    %% fuseCameras.m
    
    tic;
    
    if verbose == 2
        disp('--');
        disp(['processing slicing data for cameras ' num2str(cameras(1)) ' and ' num2str(cameras(2))]);
        disp('--');
    end;
    
    if processingMode == 0
        if globalMask(1) == 0
            averageMaskName = [sourceFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(currentSource, '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '_CHN' num2str(sChannels(2), '%.2d') '.fusionMask' inputExtension];
        else
            averageMaskName = [globalFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(globalMask(2), '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '_CHN' num2str(sChannels(2), '%.2d') '.fusionMask' inputExtension];
        end;
    else
        if globalMask(1) == 0
            averageMaskName = [sourceFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(currentSource, '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '.fusionMask' inputExtension];
        else
            averageMaskName = [globalFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(globalMask(2), '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(sChannels(1), '%.2d') '.fusionMask' inputExtension];
        end;
    end;
    averageMask = readImage(averageMaskName);
    
    if correction(2) ~= 0
        if processingMode == 0 && correction(1) ~= 0
            background = 0;
        else
            if dataType == 1 && processingMode == 2
                minIntensityName = [inputFolder '/' inputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.minIntensity.mat'];
                load(minIntensityName, 'minIntensity');
                backgroundIntensity1 = minIntensity(end);
                
                minIntensityName = [inputFolder '/' inputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.minIntensity.mat'];
                load(minIntensityName, 'minIntensity');
                backgroundIntensity2 = minIntensity(end);
            elseif dataType == 1 && processingMode == 0
                background = 0;
            else  % dataType == 0
                backgroundArray1 = primaryDataArray{1, 1}(1:subSampling:end);
                backgroundArray2 = primaryDataArray{2, 1}(1:subSampling:end);
                background = prctile(cat(2, backgroundArray1(backgroundArray1 > 0), backgroundArray2(backgroundArray2 > 0)), percentile);
                clear backgroundArray1 backgroundArray2;
            end;
        end;
        
        if correction(2) == 1
            if verbose == 2
                disp('analyzing intensity difference');
            end;
            
            minMask = min(averageMask(averageMask > 0));
            maxMask = max(averageMask(:));
            
            lowerBound = double(minMask) - double(intSizes(2));
            if lowerBound < 1
                lowerSlabSize = intSizes(2) + lowerBound - 1;
            else
                lowerSlabSize = intSizes(2);
            end;
            
            upperBound = double(maxMask) + double(intSizes(2)) - double(zSize);
            if upperBound > 0
                upperSlabSize = intSizes(2) - upperBound;
            else
                upperSlabSize = intSizes(2);
            end;
            
            finalSlabSize = lowerSlabSize + upperSlabSize + 1;
            
            dataSlice1 = zeros(xSize, zSize, finalSlabSize, 'uint16');
            dataSlice2 = zeros(xSize, zSize, finalSlabSize, 'uint16');
            for y = 1:ySize
                for x = 1:xSize
                    if averageMask(x, y) > 0
                        currentIndexArray = (averageMask(x, y) - lowerSlabSize):(averageMask(x, y) + upperSlabSize);
                        dataSlice1(x, y, :) = primaryDataArray{1, 1}(x, y, currentIndexArray);
                        dataSlice2(x, y, :) = primaryDataArray{2, 1}(x, y, currentIndexArray);
                    end;
                end;
            end;
            
            if dataType == 1 && processingMode == 2
                dataSlice1 = dataSlice1 - backgroundIntensity1;
                dataSlice2 = dataSlice2 - backgroundIntensity2;
            else % dataType == 0
                dataSlice1 = dataSlice1 - background;
                dataSlice2 = dataSlice2 - background;
            end;
            
            if processingMode == 0
                dataSlice1Name = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.fusionDataSlice' outputExtension];
                dataSlice2Name = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.fusionDataSlice' outputExtension];
            else
                dataSlice1Name = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.fusionDataSlice' outputExtension];
                dataSlice2Name = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.fusionDataSlice' outputExtension];
            end;
            
            if ~reducedIO || reducedIO == 2
                writeImage(dataSlice1, dataSlice1Name);
                writeImage(dataSlice2, dataSlice2Name);
            end;
            
            intensityList1 = double(dataSlice1(:));
            intensitySum1 = sum(intensityList1(intensityList1 > 0));
            intensityList2 = double(dataSlice2(:));
            intensitySum2 = sum(intensityList2(intensityList2 > 0));
            
            if intensitySum1 > intensitySum2
                correctionFactor = intensitySum1 / intensitySum2;
                correctionFlag = 2;
            else
                correctionFactor = intensitySum2 / intensitySum1;
                correctionFlag = 1;
            end;
            
            if dataType == 1 && processingMode == 2
                intensityCorrection = [double(backgroundIntensity1) double(backgroundIntensity2) double(intensitySum1) double(intensitySum2) double(correctionFactor) double(correctionFlag)];
            else
                intensityCorrection = [double(background) double(intensitySum1) double(intensitySum2) double(correctionFactor) double(correctionFlag)];
            end;
            intensityCorrectionName = [outputString '.intensityCorrection.mat'];
            save(intensityCorrectionName, 'intensityCorrection');
        elseif correction(2) == 2
            if correction(4) > 1
                correctionFactor = correction(4);
                correctionFlag = 2;
            else
                correctionFactor = 1 / correction(4);
                correctionFlag = 1;
            end;
        elseif correction(2) == 3
            if processingMode == 0
                correctionFactor = currentTransformations(12);
                correctionFlag = currentTransformations(13);
            else
                correctionFactor = currentTransformations(4);
                correctionFlag = currentTransformations(5);
            end;
        end;
        
        if processingMode == 0
            correctedStack1Name = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.correctedStack' outputExtension];
            correctedStack2Name = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '_CHN' num2str(tChannels(2), '%.2d') '.correctedStack' outputExtension];
        else
            correctedStack1Name = [outputHeader '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.correctedStack' outputExtension];
            correctedStack2Name = [outputHeader '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(tChannels(1), '%.2d') '.correctedStack' outputExtension];
        end;
        
        if exist(correctedStack1Name, 'file') ~= 2 || exist(correctedStack2Name, 'file') ~= 2 || enforceFlag(7)
            if verbose == 2
                disp('correcting image stacks');
            end;
            
            if dataType == 1 && processingMode == 2
                primaryDataArray{1, 1} = primaryDataArray{1, 1} - backgroundIntensity1;
                primaryDataArray{2, 1} = primaryDataArray{2, 1} - backgroundIntensity2;
            else % dataType == 0
                primaryDataArray{1, 1} = primaryDataArray{1, 1} - background;
                primaryDataArray{2, 1} = primaryDataArray{2, 1} - background;
            end;
            primaryDataArray{1, 1}(primaryDataArray{1, 1} < 0) = 0;
            primaryDataArray{2, 1}(primaryDataArray{2, 1} < 0) = 0;
            
            if correctionFlag == 1
                primaryDataArray{1, 1} = primaryDataArray{1, 1} .* correctionFactor;
            else
                primaryDataArray{2, 1} = primaryDataArray{2, 1} .* correctionFactor;
            end;
            
            if ~reducedIO
                writeImage(primaryDataArray{1, 1}, correctedStack1Name);
                writeImage(primaryDataArray{2, 1}, correctedStack2Name);
            end;
        else
            if verbose == 2
                disp('reading intensity-corrected stacks');
            end;
            
            primaryDataArray{1, 1} = readImage(correctedStack1Name);
            primaryDataArray{2, 1} = readImage(correctedStack2Name);
        end;
    else
        if verbose == 2
            disp('skipping intensity correction');
        end;
    end;
    
    if fusionType == 0
        stitchedStackName = [outputString '.stitchedStack' outputExtension];
        
        if exist(stitchedStackName, 'file') ~= 2 || enforceFlag(8)
            if verbose == 2
                disp('stitching corrected stacks');
            end;
            
            channelMask = bsxfun(@le, reshape(uint16(1:zSize), [1, 1, zSize]), averageMask);
            if frontFlag == 1
                fusedStack = bsxfun(@times, uint16(averageMask > 0), primaryDataArray{2, 1});
                fusedStack(channelMask) = primaryDataArray{1, 1}(channelMask);
            else
                fusedStack = bsxfun(@times, uint16(averageMask > 0), primaryDataArray{1, 1});
                fusedStack(channelMask) = primaryDataArray{2, 1}(channelMask);
            end;
            clear channelMask;
            
            if ~reducedIO
                writeImage(fusedStack, stitchedStackName);
            end;
        else
            if verbose == 2
                disp('reading stitched stack');
            end;
            
            fusedStack = readImage(stitchedStackName);
        end;
    end;
    
    fusedStackName = [outputString '.fusedStack' outputExtension];
    
    if exist(fusedStackName, 'file') ~= 2 || enforceFlag(9)
        if fusionType == 0
            if verbose == 2
                disp('blending corrected stacks');
            end;
            
            weighting1Front = reshape((1:((0.5 - 1) / (blending(2) - 1)):0.5), [1 1 blending(2)]); % dominant fraction
            weighting2Front = reshape((0:((0.5 - 0) / (blending(2) - 1)):0.5), [1 1 blending(2)]); % fading fraction
            
            weighting1Back = reshape((0.5:((1 - 0.5) / (blending(2) - 1)):1), [1 1 blending(2)]); % dominant fraction
            weighting2Back = reshape((0.5:((0 - 0.5) / (blending(2) - 1)):0), [1 1 blending(2)]); % fading fraction
            
            if frontFlag == 1
                for y = 1:ySize
                    for x = 1:xSize
                        if averageMask(x, y) > 0
                            array1Front = double(primaryDataArray{1, 1}(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)));
                            array2Front = double(primaryDataArray{2, 1}(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)));
                            array1Back = double(primaryDataArray{1, 1}(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))));
                            array2Back = double(primaryDataArray{2, 1}(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))));
                            fusedStack(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)) = uint16(array1Front .* weighting1Front + array2Front .* weighting2Front);
                            fusedStack(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))) = uint16(array1Back .* weighting2Back + array2Back .* weighting1Back);
                        end;
                    end;
                end;
            else
                for y = 1:ySize
                    for x = 1:xSize
                        if averageMask(x, y) > 0
                            array1Front = double(primaryDataArray{1, 1}(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)));
                            array2Front = double(primaryDataArray{2, 1}(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)));
                            array1Back = double(primaryDataArray{1, 1}(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))));
                            array2Back = double(primaryDataArray{2, 1}(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))));
                            fusedStack(x, y, (averageMask(x, y) - blending(2) + 1):averageMask(x, y)) = uint16(array1Front .* weighting2Front + array2Front .* weighting1Front);
                            fusedStack(x, y, (averageMask(x, y) + 1):(averageMask(x, y) + blending(2))) = uint16(array1Back .* weighting1Back + array2Back .* weighting2Back);
                        end;
                    end;
                end;
            end;
        elseif fusionType == 1
            if verbose == 2
                disp('performing wavelet fusion of corrected stacks');
            end;
            
            fusedStack = zeros(xSize, ySize, zSize, 'uint16');
            for z = 1:zSize
                fusedStack(:, :, z) = uint16(wfusimg(primaryDataArray{1, 1}(:, :, z), primaryDataArray{2, 1}(:, :, z), 'db4', 5, 'mean', 'max'));
            end;
        elseif fusionType == 2
            if verbose == 2
                disp('performing averaging of corrected stacks');
            end;
            
            fusedStack = zeros(xSize, ySize, zSize, 'uint16');
            for z = 1:zSize
                fusedStack(:, :, z) = (primaryDataArray{1, 1}(:, :, z) + primaryDataArray{2, 1}(:, :, z)) ./ 2;
            end;
        end;
        
        writeImage(fusedStack, fusedStackName);
        
        fusedStackProjXYName = [outputString '.fusedStack_xyProjection' outputExtension];
        fusedStackProjXZName = [outputString '.fusedStack_xzProjection' outputExtension];
        fusedStackProjYZName = [outputString '.fusedStack_yzProjection' outputExtension];
        
        writeImage(max(fusedStack, [], 3), fusedStackProjXYName);
        writeImage(squeeze(max(fusedStack, [], 2)), fusedStackProjXZName);
        writeImage(squeeze(max(fusedStack, [], 1)), fusedStackProjYZName);
    end;
    
    clear fusedStack;
    
    if verbose == 2
        disp(' ');
    end;
    
    time = toc;
    if verbose == 1
        disp(['***     camera fusion completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
        disp(' ');
    elseif verbose == 2
        disp(['*** camera fusion completed in ' num2str(floor(floor(time / 60) / 60), '%.2d') ':' num2str(mod(floor(time / 60), 60), '%.2d') ':' num2str(round(mod(time, 60)), '%.2d')]);
        disp(' ');
    end;
end;

jobCompleted = 1;
save([outputString '_jobCompleted.txt'], '-ascii', 'jobCompleted');

end