timepoints = 0:278;

fusionAxis = 3; % 2: y-axis (SPM1), 3: z-axis (SPM0)
headers = {...
    'SPM00_TM';...
    'SPM00_TM'};

% fusionAxis = 2; % 2: y-axis (SPM1), 3: z-axis (SPM0)
% headers = {...
%     'SPM01_TM';...
%     'SPM01_TM'};

footers = {...
    '_CM00_CHN00.affine.trsf.cropped.padded.cropped.klb';...
    '_CM01_CHN00.affine.trsf.cropped.padded.cropped.klb'};

splitting   = 10;
kernelSize  = 5;
kernelSigma = 2;

background  = 100;
threshold   = 0.5;

blending    = 20;

projectionFolder = 'ProjectionsBlending';

poolWorkers = 6;

%% main loop

if exist(projectionFolder, 'dir') ~= 7
    mkdir(projectionFolder);
end;

if exist([projectionFolder '\Frames'], 'dir') ~= 7
    mkdir([projectionFolder '\Frames']);
end;

kernelSizeArray  = [1 1 1] .* kernelSize;
kernelSigmaArray = [1 1 1] .* kernelSigma;

disp(' ');
if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);
disp(' ');

parfor tIndex = 1:numel(timepoints)
    t = timepoints(tIndex);
    
    if exist([projectionFolder '\Frames\' headers{1} num2str(t, '%.6d') '.yz.klb'], 'file') == 2
        disp(['skipping time point ' num2str(t)]);
    else
        disp(['processing time point ' num2str(t)]);
        
        stack1 = readImage(['TM' num2str(t, '%.6d') '\' headers{1} num2str(t, '%.6d') footers{1}]);
        stack2 = readImage(['TM' num2str(t, '%.6d') '\' headers{2} num2str(t, '%.6d') footers{2}]);
        
        if fusionAxis == 3
            stack1 = permute(stack1, [1 3 2]);
            stack2 = permute(stack2, [1 3 2]);
        end;
        
        sumStack = (stack1 + stack2) / 2;
        
        [xSize, ySize, zSize] = size(sumStack);
        
        if splitting > 1
            gaussStack = zeros(xSize, ySize, zSize, 'uint16');
            
            splittingMargin = 2 * kernelSize;
            
            for i = 1:splitting
                xSlabStart = max(1, round((i - 1) * xSize / splitting + 1 - splittingMargin));
                xSlabStop = min(xSize, round(i * xSize / splitting + splittingMargin));
                convolvedSlab = uint16(imgaussianAnisotropy(double(sumStack(xSlabStart:xSlabStop, :, :)), kernelSigmaArray, kernelSizeArray));
                if i == 1
                    gaussStack(1:(xSlabStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                elseif i == splitting
                    gaussStack((xSlabStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                else % i > 1 && i < splitting
                    gaussStack((xSlabStart + splittingMargin):(xSlabStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                end;
            end;
        else
            gaussStack = uint16(imgaussianAnisotropy(double(sumStack), kernelSigmaArray, kernelSizeArray));
        end;
        
        intensityStatistics = zeros(splitting, 2);
        for i = 1:splitting
            xSlabStart = round((i - 1) * xSize / splitting + 1);
            xSlabStop  = round(i * xSize / splitting);
            temporaryArray = gaussStack(xSlabStart:xSlabStop, :, :);
            temporaryArray = temporaryArray(temporaryArray > background);
            intensityStatistics(i, 1) = sum(temporaryArray(:));
            intensityStatistics(i, 2) = size(temporaryArray, 1);
        end;
        meanIntensity = sum(intensityStatistics(:, 1)) / sum(intensityStatistics(:, 2));
        
        level = background + (meanIntensity - background) * threshold;
        
        mask = gaussStack > level;
        
        transitionMask = zeros(xSize, zSize, 'uint16');
        for i = 1:splitting
            xSlabStart = round((i - 1) * xSize / splitting + 1);
            xSlabStop  = round(i * xSize / splitting);
            coordinateMask = repmat(1:ySize, [(xSlabStop - xSlabStart + 1) 1 size(mask, 3)]);
            coordinateMask(mask(xSlabStart:xSlabStop, :, :) == 0) = NaN;
            transitionMask(xSlabStart:xSlabStop, :) = uint16(round(squeeze(nanmean(coordinateMask, 2))));
        end;
        
        transitionMask(transitionMask > (ySize - blending)) = round(ySize / 2);
        transitionMask(transitionMask < blending) = round(ySize / 2);
        transitionMask(transitionMask == 0) = round(ySize / 2);
        
        writeImage(transitionMask, ['TM' num2str(t, '%.6d') '\' headers{1} num2str(t, '%.6d') '.transitionMask.klb']);
        
        stitchingMask = bsxfun(@le, reshape(uint16(1:ySize), [1, ySize, 1]), reshape(transitionMask, [xSize, 1, zSize]));
        fusedStack = stack2;
        fusedStack(stitchingMask) = stack1(stitchingMask);
        
        weighting1Left = 1:((0.5 - 1) / (blending(1) - 1)):0.5; % dominant fraction
        weighting2Left = 0:((0.5 - 0) / (blending(1) - 1)):0.5; % fading fraction
        
        weighting1Right = 0.5:((1 - 0.5) / (blending(1) - 1)):1; % dominant fraction
        weighting2Right = 0.5:((0 - 0.5) / (blending(1) - 1)):0; % fading fraction
        
        for z = 1:zSize
            for x = 1:xSize
                if transitionMask(x, z) > 0
                    array1Left = double(stack1(x, (transitionMask(x, z) - blending(1) + 1):transitionMask(x, z), z));
                    array2Left = double(stack2(x, (transitionMask(x, z) - blending(1) + 1):transitionMask(x, z), z));
                    array1Right = double(stack1(x, (transitionMask(x, z) + 1):(transitionMask(x, z) + blending(1)), z));
                    array2Right = double(stack2(x, (transitionMask(x, z) + 1):(transitionMask(x, z) + blending(1)), z));
                    fusedStack(x, (transitionMask(x, z) - blending(1) + 1):transitionMask(x, z), z) = uint16(array1Left .* weighting1Left + array2Left .* weighting2Left);
                    fusedStack(x, (transitionMask(x, z) + 1):(transitionMask(x, z) + blending(1)), z) = uint16(array1Right .* weighting2Right + array2Right .* weighting1Right);
                end;
            end;
        end;
        
        fusedStack = permute(fusedStack, [1 3 2]);
        
        writeImage(fusedStack, ['TM' num2str(t, '%.6d') '\' headers{1} num2str(t, '%.6d') '.fusedStack.klb']);
        
        writeImage(max(fusedStack, [], 3), [projectionFolder '\Frames\' headers{1} num2str(t, '%.6d') '.xy.klb']);
        writeImage(squeeze(max(fusedStack, [], 2)), [projectionFolder '\Frames\' headers{1} num2str(t, '%.6d') '.xz.klb']);
        writeImage(squeeze(max(fusedStack, [], 1)), [projectionFolder '\Frames\' headers{1} num2str(t, '%.6d') '.yz.klb']);
    end;
end;

disp(' ');
if matlabpool('size') > 0
    matlabpool('close');
end;
disp(' ');