inputFolder      = 'X:' filesep 'SV1' filesep '14-05-21' filesep 'Mmu_E1_CAGTAG1.corrected' filesep 'Results' filesep 'TimeFused.Corrected';
inputHeader      = 'Mmu_E1_CAGTAG1.TM';
inputFooter      = '_timeFused_blending';
inputFileHeader  = 'SPM00_TM';
inputFileMiddle  = '_CM00_CM01_CHN';
inputFileFooter  = '.fusedStack.corrected.shifted.klb';

isotropicFolder  = 'X:' filesep 'SV1' filesep '14-05-21' filesep 'Mmu_E1_CAGTAG1.corrected' filesep 'Results' filesep 'TimeFused.Corrected.Interpolated';
masksFolder      = 'X:' filesep 'SV1' filesep '14-05-21' filesep 'Mmu_E1_CAGTAG1.corrected' filesep 'Results' filesep 'Peeling' filesep 'Masks.T9.ThetaAP';
shellsFolder     = 'X:' filesep 'SV1' filesep '14-05-21' filesep 'Mmu_E1_CAGTAG1.corrected' filesep 'Results' filesep 'Peeling' filesep 'Shells.T9.ThetaAP';

timePoints       = 0:531;
masterChannels   = 1;
childChannels    = 0;

embryoCenters    = [...
      0 1450  945 477;...
     50 1425  925 479;...
    100 1500  925 479;...
    150 1650  970 472;...
    200 1800 1000 472;...
    531 1800 1000 472];

axialScaling     = 2.031/(6.5/16);
splitFactor      = 20;
showProjections  = 0;

kernelSize       = 5;
kernelSigma      = 2;
minIntensity     = 0;
thresholdFactor  = 0.9;
fraction         = 10^-5;

coordinateSystem = 1;       % see definitions below
rCoarseStepSize  = 1;       % in pixels/voxels
rFineStepSize    = 0.5;     % in pixels/voxels
thetaSteps       = 3000;
phiSteps         = 3000;

medianType       = 0;       % 0: medfilt2, 1:mediannan
fillGaps         = [1 150]; % slot 1: yes/no, slot 2: size of median filter for gap replacement

removeBumps      = [2 2];   % slot 1: 0: no bump removal, 1: conservative bump removal, 2: agressive bump removal (may cause some geometrical artifacts)
                            % slot 2: number of smoothing iterations

shellThicknesses = [...     % in pixels/voxels
      0 10;...
    100 10;...
    200 30;...
    531 30];

save3DShells     = 0;
saveOrthoMaps    = 0;

% Coordinate system #0:
% x =  r sin_theta cos_phi
% y = -r cos_theta
% z =  r sin_theta sin_phi

% Coordinate system #1:
% x =  r sin_theta cos_phi
% y = -r sin_theta sin_phi
% z = -r cos_theta

if exist(isotropicFolder, 'dir') ~= 7
    mkdir(isotropicFolder);
end;

if exist(masksFolder, 'dir') ~= 7
    mkdir(masksFolder);
end;

if exist(shellsFolder, 'dir') ~= 7
    mkdir(shellsFolder);
end;

thetaRange = 0:(pi/thetaSteps):pi;
phiRange   = 0:(2*pi/phiSteps):2*pi;

kernelSizeArray  = [kernelSize kernelSize kernelSize];
kernelSigmaArray = [kernelSigma kernelSigma kernelSigma];

embryoCenterArray = zeros(numel(embryoCenters(1, 1):embryoCenters(end, 1)), 4);
embryoCenterArray(:, 1) = embryoCenters(1, 1):embryoCenters(end, 1);

for i = 1:(size(embryoCenters, 1) - 1)
    currentArray = zeros(numel(embryoCenters(i, 1):embryoCenters(i + 1, 1)), 4);
    currentArray(:, 1) = embryoCenters(i, 1):embryoCenters(i + 1, 1);
    if embryoCenters(i, 2) == embryoCenters(i + 1, 2)
         currentArray(:, 2) = embryoCenters(i, 2);
    else
        currentArray(:, 2) = embryoCenters(i, 2):(embryoCenters(i + 1, 2) - embryoCenters(i, 2))/(numel(embryoCenters(i, 1):embryoCenters(i + 1, 1)) - 1):embryoCenters(i + 1, 2);
    end;
    if embryoCenters(i, 3) == embryoCenters(i + 1, 3)
         currentArray(:, 3) = embryoCenters(i, 3);
    else
        currentArray(:, 3) = embryoCenters(i, 3):(embryoCenters(i + 1, 3) - embryoCenters(i, 3))/(numel(embryoCenters(i, 1):embryoCenters(i + 1, 1)) - 1):embryoCenters(i + 1, 3);
    end;
    if embryoCenters(i, 4) == embryoCenters(i + 1, 4)
         currentArray(:, 4) = embryoCenters(i, 4);
    else
        currentArray(:, 4) = embryoCenters(i, 4):(embryoCenters(i + 1, 4) - embryoCenters(i, 4))/(numel(embryoCenters(i, 1):embryoCenters(i + 1, 1)) - 1):embryoCenters(i + 1, 4);
    end;
    for j = 1:size(currentArray, 1)
        embryoCenterArray(find(embryoCenterArray(:, 1) == currentArray(j, 1), 1), :) = currentArray(j, :);
    end;
end;

embryoCenterArray(:, 4) = embryoCenterArray(:, 4) * axialScaling;
embryoCenterArray(:, 2:4) = round(embryoCenterArray(:, 2:4));

shellThicknessArray = zeros(numel(shellThicknesses(1, 1):shellThicknesses(end, 1)), 2);
shellThicknessArray(:, 1) = shellThicknesses(1, 1):shellThicknesses(end, 1);

for i = 1:(size(shellThicknesses, 1) - 1)
    currentArray = zeros(numel(shellThicknesses(i, 1):shellThicknesses(i + 1, 1)), 2);
    currentArray(:, 1) = shellThicknesses(i, 1):shellThicknesses(i + 1, 1);
    if shellThicknesses(i, 2) == shellThicknesses(i + 1, 2)
         currentArray(:, 2) = shellThicknesses(i, 2);
    else
        currentArray(:, 2) = shellThicknesses(i, 2):(shellThicknesses(i + 1, 2) - shellThicknesses(i, 2))/(numel(shellThicknesses(i, 1):shellThicknesses(i + 1, 1)) - 1):shellThicknesses(i + 1, 2);
    end;
    for j = 1:size(currentArray, 1)
        shellThicknessArray(find(shellThicknessArray(:, 1) == currentArray(j, 1), 1), :) = currentArray(j, :);
    end;
end;

shellThicknessArray(:, 2) = round(shellThicknessArray(:, 2));

for tIndex = 1:numel(timePoints)
    t = timePoints(tIndex);
    
    embryoCenter = embryoCenterArray(find(embryoCenterArray(:, 1) == t, 1), 2:4);
    rShellThickness = shellThicknessArray(find(shellThicknessArray(:, 1) == t, 1), 2);
    
    if fillGaps(1)
        medianString = ['.filledGaps_' num2str(fillGaps(2))];
    else
        medianString = '';
    end;
    
    if removeBumps(1) > 0
        smoothString = ['.smoothed_' num2str(removeBumps(1)) '_i' num2str(removeBumps(2))];
    else
        smoothString = '';
    end;
    
    lastMasterStackName = [inputFileHeader num2str(t, '%.6d') inputFileMiddle num2str(masterChannels(end), '%.2d') inputFileFooter];
    lastMasterShellProjectionName = [lastMasterStackName(1:(end-3)) 'shellProjection_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
    
    lastChildStackName = [inputFileHeader num2str(t, '%.6d') inputFileMiddle num2str(childChannels(end), '%.2d') inputFileFooter];
    lastChildShellProjectionName = [lastChildStackName(1:(end-3)) 'shellProjection_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
    
    if exist([shellsFolder filesep '' lastMasterShellProjectionName], 'file') == 2 && exist([shellsFolder filesep '' lastChildShellProjectionName], 'file') == 2
        disp(['Last master/child shell projections detected for time point ' num2str(t)]);
    else    
        for cIndex = 1:numel(masterChannels)
            
            %% create/load interpolated stack(s) for master channel(s)
            
            c = masterChannels(cIndex);
            
            inputStackName = [inputFileHeader num2str(t, '%.6d') inputFileMiddle num2str(c, '%.2d') inputFileFooter];
            outputInterpolatedStackName = [inputStackName(1:(end-3)) 'interpolated.klb'];
            
            inputFullPath = [inputFolder filesep '' inputHeader num2str(t, '%.6d') inputFooter filesep '' inputStackName];
            outputFullPath = [isotropicFolder filesep '' outputInterpolatedStackName];
            
            if exist(outputFullPath, 'file') == 2
                disp(['Existing interpolated stack detected for time point ' num2str(t) ', channel ' num2str(c)]);
                
                interpolatedStack = readImage(outputFullPath);
            else
                stack = readImage(inputFullPath);
                
                if showProjections
                    figure; imagesc(max(stack, [], 3));
                    figure; imagesc(squeeze(max(stack, [], 2)));
                end;
                
                if splitFactor == 1
                    xArray = 1:size(stack, 1);
                    yArray = 1:size(stack, 2);
                    zArray = 1:size(stack, 3);
                    zArrayInterpolated = 1:(1/axialScaling):size(stack, 3);
                    
                    xArrayLength = numel(xArray);
                    yArrayLength = numel(yArray);
                    zArrayLength = numel(zArray);
                    zArrayInterpolatedLength = numel(zArrayInterpolated);
                    
                    interpolatedStack = uint16(interp3(...
                        repmat(yArray, xArrayLength, 1, zArrayLength), ...
                        repmat(xArray', 1, yArrayLength, zArrayLength), ...
                        repmat(permute(zArray, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                        single(stack), ...
                        repmat(yArray, xArrayLength, 1, zArrayInterpolatedLength), ...
                        repmat(xArray', 1, yArrayLength, zArrayInterpolatedLength), ...
                        repmat(permute(zArrayInterpolated, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                        'cubic'));
                else
                    zArray = 1:size(stack, 3);
                    zArrayInterpolated = 1:(1/axialScaling):size(stack, 3);
                    
                    zArrayLength = numel(zArray);
                    zArrayInterpolatedLength = numel(zArrayInterpolated);
                    
                    interpolatedStack = zeros(size(stack, 1), size(stack, 2), zArrayInterpolatedLength, 'uint16');
                    
                    for x = 1:splitFactor
                        xStart = (x - 1) * ceil(size(stack, 1) / splitFactor) + 1;
                        if x == splitFactor
                            xStop = size(stack, 1);
                        else
                            xStop = xStart + ceil(size(stack, 1) / splitFactor) - 1;
                        end;
                        
                        xArray = 1:(xStop - xStart + 1);
                        xArrayLength = numel(xArray);
                        
                        for y = 1:splitFactor
                            disp(['Processing slab (' num2str(x) ', ' num2str(y) ') for time point ' num2str(t) ', channel ' num2str(c)]);
                            
                            yStart = (y - 1) * ceil(size(stack, 2) / splitFactor) + 1;
                            if y == splitFactor
                                yStop = size(stack, 2);
                            else
                                yStop = yStart + ceil(size(stack, 2) / splitFactor) - 1;
                            end;
                            
                            yArray = 1:(yStop - yStart + 1);
                            yArrayLength = numel(yArray);
                            
                            interpolatedStack(xStart:xStop, yStart:yStop, :) = uint16(interp3(...
                                repmat(yArray, xArrayLength, 1, zArrayLength), ...
                                repmat(xArray', 1, yArrayLength, zArrayLength), ...
                                repmat(permute(zArray, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                single(stack(xStart:xStop, yStart:yStop, :)), ...
                                repmat(yArray, xArrayLength, 1, zArrayInterpolatedLength), ...
                                repmat(xArray', 1, yArrayLength, zArrayInterpolatedLength), ...
                                repmat(permute(zArrayInterpolated, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                'cubic'));
                        end;
                    end;
                end;
                
                writeImage(interpolatedStack, outputFullPath);
                
                if showProjections
                    figure; imagesc(max(interpolatedStack, [], 3));
                    figure; imagesc(squeeze(max(interpolatedStack, [], 2)));
                end;
            end;
            
            %% create/load adaptively thresholded mask(s) for master channel(s)
            
            [xSize, ySize, zSize] = size(interpolatedStack);
            
            outputMaskName = [inputStackName(1:(end-3)) 'mask.klb'];
            outputFullPath = [masksFolder filesep '' outputMaskName];
            
            if exist(outputFullPath, 'file') == 2
                disp(['Existing mask detected for time point ' num2str(t) ', channel ' num2str(c)]);
                
                mask = readImage(outputFullPath) > 0;
            else
                if splitFactor > 1
                    gaussStack = zeros(xSize, ySize, zSize, 'uint16');
                    
                    splittingMargin = 2 * kernelSize;
                    
                    for i = 1:splitFactor
                        disp(['Convolving slab ' num2str(i) ' for time point ' num2str(t) ', channel ' num2str(c)]);
                        
                        xSlabStart = max(1, round((i - 1) * xSize / splitFactor + 1 - splittingMargin));
                        xSlabStop = min(xSize, round(i * xSize / splitFactor + splittingMargin));
                        convolvedSlab = imgaussianAnisotropy(interpolatedStack(xSlabStart:xSlabStop, :, :), kernelSigmaArray, kernelSizeArray);
                        if i == 1
                            gaussStack(1:(xSlabStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                        elseif i == splitFactor
                            gaussStack((xSlabStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                        else % i > 1 && i < splitFactor
                            gaussStack((xSlabStart + splittingMargin):(xSlabStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                        end;
                        clear convolvedSlab;
                    end;
                else
                    gaussStack = imgaussianAnisotropy(stack, kernelSigmaArray, kernelSizeArray);
                end;
                
                meanIntensity = mean(gaussStack(gaussStack > 0));
                level = minIntensity + (meanIntensity - minIntensity) * thresholdFactor;
                
                mask = gaussStack > level;
                if fraction > 0
                    threshold = round(xSize * ySize * zSize * fraction);
                    mask = bwareaopen(mask, threshold);
                end;
                
                writeImage(uint8(mask), outputFullPath);
            end;
            
            %% create/load shell map(s) for master channel(s)
            
            phiSize = numel(phiRange);
            thetaSize = numel(thetaRange);
            
            outputShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) '.mat'];
            outputFullPath = [shellsFolder filesep '' outputShellMapName];
            
            if exist(outputFullPath, 'file') == 2
                disp(['Existing shell map detected for time point ' num2str(t) ', channel ' num2str(c)]);
                
                load(outputFullPath);
            else
                disp(['Creating shell map for time point ' num2str(t) ', channel ' num2str(c)]);
                
                maxDistanceX = max(embryoCenter(1), xSize - embryoCenter(1));
                maxDistanceY = max(embryoCenter(2), ySize - embryoCenter(2));
                maxDistanceZ = max(embryoCenter(3), zSize - embryoCenter(3));
                maxR = ceil(sqrt(maxDistanceX^2 + maxDistanceY^2 + maxDistanceZ^2));
                
                rRangeCoarse = 0:rCoarseStepSize:maxR;
                
                shellMap = NaN(phiSize, thetaSize);
                
                for i = 1:numel(phiRange)
                    for j = 1:numel(thetaRange)
                        switch coordinateSystem
                            case 0
                                xArray = round(embryoCenter(1) + rRangeCoarse .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                yArray = round(embryoCenter(2) - rRangeCoarse .* cos(thetaRange(j)));
                                zArray = round(embryoCenter(3) + rRangeCoarse .* (sin(thetaRange(j)) * sin(phiRange(i))));
                            case 1
                                xArray = round(embryoCenter(1) + rRangeCoarse .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                yArray = round(embryoCenter(2) - rRangeCoarse .* (sin(thetaRange(j)) * sin(phiRange(i))));
                                zArray = round(embryoCenter(3) - rRangeCoarse .* cos(thetaRange(j)));
                        end;
                        
                        outOfBoundsFlag = xArray < 1 | xArray > xSize | yArray < 1 | yArray > ySize | zArray < 1 | zArray > zSize;
                        
                        linearIndices = sub2ind(size(mask), xArray(~outOfBoundsFlag), yArray(~outOfBoundsFlag), zArray(~outOfBoundsFlag));
                        mostRemoteVoxel = find(mask(linearIndices), 1, 'last');
                        if ~isempty(mostRemoteVoxel)
                            shellMap(i, j) = rRangeCoarse(mostRemoteVoxel);
                        end;
                    end;
                end;
                
                % figure;
                % h = imagesc(shellMap);
                % xlabel('theta');
                % ylabel('phi');
                % set(gca, 'YTick', [1 phiSize]);
                % set(gca, 'YTickLabel', [phiRange(1) phiRange(phiSize)]);
                % set(gca, 'XTick', [1 thetaSize]);
                % set(gca, 'XTickLabel', [thetaRange(1) thetaRange(thetaSize)]);
                % colorbar;
                
                save(outputFullPath, 'shellMap');
                
                outputShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellMapName];
                
                writeImage(shellMap, outputFullPath);
            end;
            
            %% (optional) replace NaNs in shell map with median of local neighborhood
            
            if fillGaps(1)
                medianString = ['.filledGaps_' num2str(fillGaps(2))];
                
                outputPatchedShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) medianString '.mat'];
                outputFullPath = [shellsFolder filesep '' outputPatchedShellMapName];
                
                if exist(outputFullPath, 'file') == 2
                    disp(['Existing gap-filled shell map detected for time point ' num2str(t) ', channel ' num2str(c)]);
                    
                    load(outputFullPath);
                else
                    disp(['Creating gap-filled shell map for time point ' num2str(t) ', channel ' num2str(c)]);
                    
                    filteredShellMap = shellMap;
                    filteredShellMap(isnan(filteredShellMap)) = 0;
                    
                    if medianType == 0
                        filteredShellMap = medfilt2(filteredShellMap, [fillGaps(2) fillGaps(2)]);
                    else
                        filteredShellMap = mediannan(filteredShellMap, [fillGaps(2) fillGaps(2)]);
                    end;
                    
                    replacementFlag = isnan(shellMap) & filteredShellMap ~= 0;
                    
                    shellMap(replacementFlag) = filteredShellMap(replacementFlag);
                    
                    save(outputFullPath, 'shellMap');
                    
                    outputPatchedShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) medianString '.klb'];
                    outputFullPath = [shellsFolder filesep '' outputPatchedShellMapName];
                    
                    writeImage(shellMap, outputFullPath);
                end;
            else
                medianString = '';
            end;
            
            %% (optional) remove positive and negative bumps in shell map
            
            if removeBumps(1) > 0
                for smoothingIteration = 1:removeBumps(2)
                    smoothString = ['.smoothed_' num2str(removeBumps(1)) '_i' num2str(smoothingIteration)];
                    
                    outputSmoothShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) medianString smoothString '.mat'];
                    outputFullPath = [shellsFolder filesep '' outputSmoothShellMapName];
                    
                    if exist(outputFullPath, 'file') == 2
                        disp(['Existing smoothed shell map (i' num2str(smoothingIteration) ') detected for time point ' num2str(t) ', channel ' num2str(c)]);
                        
                        load(outputFullPath);
                    else
                        disp(['Creating smoothed shell map (i' num2str(smoothingIteration) ') for time point ' num2str(t) ', channel ' num2str(c)]);
                        
                        % Step 1: perform edge detection on shell map
                        shellMapForEdgeDetection = uint16(shellMap);
                        shellMapForEdgeDetection(isnan(shellMap)) = 0;
                        edgeMap = edge(shellMapForEdgeDetection, 'canny');
                        % figure; imagesc(edgeMap)
                        
                        structuralElement = strel('disk', 2);
                        erodedBinaryMap = imerode(shellMapForEdgeDetection > 0, structuralElement);
                        erodedEdgeMap = edgeMap;
                        erodedEdgeMap(~erodedBinaryMap) = 0;
                        % figure; imagesc(erodedEdgeMap)
                        
                        connectedComponents = bwconncomp(erodedEdgeMap);
                        % figure; imagesc(imfill(labelMatrix(connectedComponents), 'holes'));
                        % figure; imagesc(imfill(labelMatrix(connectedComponents) > 0, 'holes'));
                        
                        % Step 2: close connected components by iterative dilations
                        blankMap = false(connectedComponents.ImageSize);
                        regionStatistics = regionprops(connectedComponents, 'ConvexImage', 'EulerNumber');
                        for i = find([regionStatistics.EulerNumber] > 0)
                            distanceImage = bwdist(~regionStatistics(i).ConvexImage);
                            maxDistance = ceil(max(distanceImage(:)));
                            currentSlice = blankMap;
                            currentSlice(connectedComponents.PixelIdxList{i}) = true;
                            if isinf(maxDistance)
                                continue;
                            end;
                            for dilationSize = 2:maxDistance
                                newSlice = imdilate(currentSlice, ones(dilationSize));
                                regionStatisticsNew = regionprops(newSlice, 'EulerNumber');
                                if regionStatisticsNew.EulerNumber <= 0
                                    newSlice = imerode(imfill(newSlice, 'holes'), ones(dilationSize));
                                    connectedComponents.PixelIdxList{i} = find(newSlice);
                                end;
                            end;
                        end;
                        % figure; imagesc(imfill(labelmatrix(connectedComponents), 'holes'));
                        % figure; imagesc(imfill(labelmatrix(connectedComponents) > 0, 'holes'));
                        
                        % Step 3: if dilation filled everything in step 2, replace respective object with convex hull
                        if removeBumps(1) > 1
                            regionStatistics = regionprops(connectedComponents, 'ConvexImage', 'EulerNumber', 'BoundingBox');
                            for i = find([regionStatistics.EulerNumber] > 0)
                                maxDistance = ceil(max(distanceImage(:)));
                                currentSlice = blankMap;
                                currentSlice(connectedComponents.PixelIdxList{i}) = true;
                                distanceImage = bwdist(~currentSlice);
                                if ~any(distanceImage(:) > 1)
                                    newSlice = currentSlice;
                                    boundingBox = ceil(regionStatistics(i).BoundingBox);
                                    newSlice(...
                                        (1:boundingBox(4)) + boundingBox(2) - 1,...
                                        (1:boundingBox(3)) + boundingBox(1) - 1) = regionStatistics(i).ConvexImage;
                                    connectedComponents.PixelIdxList{i} = find(newSlice);
                                end;
                            end;
                        end;
                        
                        filledEdgeMap = imfill(labelmatrix(connectedComponents), 'holes');
                        % figure; imagesc(filledEdgeMap);
                        
                        % Step 4: replace indices of enclosed objects with the dominant surrounding object index
                        indicesOfHoles = find(arrayfun(@(i)mode(double(filledEdgeMap(bwmorph(filledEdgeMap == i, 'dilate', 1) & ~(filledEdgeMap == i)))), 1:connectedComponents.NumObjects));
                        for i = indicesOfHoles
                            surroundingsBinary = bwmorph(filledEdgeMap == i, 'dilate', 1) & ~(filledEdgeMap == i);
                            surroundings = filledEdgeMap(surroundingsBinary);
                            filledEdgeMap(filledEdgeMap == i) = mode(surroundings(surroundings > 0));
                        end;
                        % figure; imagesc(filledEdgeMap);
                        
                        % Step 5: remove small structures attached to blob-shaped primary objects
                        cleanEdgeMap = bwmorph(imopen(filledEdgeMap, structuralElement), 'dilate', 1);
                        % figure; imagesc(cleanEdgeMap);
                        
                        % Step 6: use processed edge map to replace flagged areas in original shell map
                        [A, B] = meshgrid(1:size(shellMap, 2), 1:size(shellMap, 1));
                        foregroundSlots = find(shellMap > 0 & ~(cleanEdgeMap > 0));
                        transitionSlots = find(cleanEdgeMap > 0);
                        f = scatteredInterpolant(A(foregroundSlots), B(foregroundSlots), double(shellMap(foregroundSlots)), 'natural');
                        transitionValues = f(A(transitionSlots), B(transitionSlots));
                        shellMap(transitionSlots) = transitionValues;
                        % figure; imagesc(shellMap);
                        
                        shellMap = round(shellMap);
                        
                        save(outputFullPath, 'shellMap');
                        
                        outputSmoothShellMapName = [inputStackName(1:(end-3)) 'shellMap_' num2str(phiSteps) '_' num2str(thetaSteps) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputSmoothShellMapName];
                        
                        writeImage(shellMap, outputFullPath);
                    end;
                end;
            else
                smoothString = '';
            end;
            
            %% use shell map(s) to create shell projection(s)/stack(s) for master channel(s)
            
            disp(['Creating shell projection for time point ' num2str(t) ', channel ' num2str(c)]);
            
            minMap = shellMap - rShellThickness;
            maxMap = shellMap;
            
            maxMap(isnan(minMap) | minMap < 1) = NaN;
            minMap(isnan(minMap) | minMap < 1) = NaN;
            
            shellProjection = NaN(phiSize, thetaSize);
            
            if save3DShells || saveOrthoMaps
                shellStack = zeros(size(interpolatedStack, 1), size(interpolatedStack, 2), size(interpolatedStack, 3), 'uint16');
            end;
            
            for i = 1:numel(phiRange)
                for j = 1:numel(thetaRange)
                    if ~isnan(minMap(i, j))
                        rRangeLocal = minMap(i, j):rFineStepSize:maxMap(i, j);
                        
                        switch coordinateSystem
                            case 0
                                xArray = round(embryoCenter(1) + rRangeLocal .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                yArray = round(embryoCenter(2) - rRangeLocal .* cos(thetaRange(j)));
                                zArray = round(embryoCenter(3) + rRangeLocal .* (sin(thetaRange(j)) * sin(phiRange(i))));
                            case 1
                                xArray = round(embryoCenter(1) + rRangeLocal .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                yArray = round(embryoCenter(2) - rRangeLocal .* (sin(thetaRange(j)) * sin(phiRange(i))));
                                zArray = round(embryoCenter(3) - rRangeLocal .* cos(thetaRange(j)));
                        end;
                        
                        outOfBoundsFlag = xArray < 1 | xArray > xSize | yArray < 1 | yArray > ySize | zArray < 1 | zArray > zSize;
                        
                        xRegion = xArray(~outOfBoundsFlag);
                        yRegion = yArray(~outOfBoundsFlag);
                        zRegion = zArray(~outOfBoundsFlag);
                        
                        if ~isempty(xRegion)
                            linearIndices = sub2ind(size(mask), xRegion, yRegion, zRegion);
                            shellProjection(i, j) = max(interpolatedStack(linearIndices));
                            
                            if save3DShells || saveOrthoMaps
                                shellStack(linearIndices) = interpolatedStack(linearIndices);
                            end;
                        end;
                    end;
                end;
            end;
            
            % figure;
            % imagesc(shellProjection);
            % xlabel('theta');
            % ylabel('phi');
            % set(gca, 'YTick', [1 phiSize]);
            % set(gca, 'YTickLabel', [phiRange(1) phiRange(phiSize)]);
            % set(gca, 'XTick', [1 thetaSize]);
            % set(gca, 'XTickLabel', [thetaRange(1) thetaRange(thetaSize)]);
            % colorbar;
            
            outputShellProjectionName = [inputStackName(1:(end-3)) 'shellProjection_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
            outputFullPath = [shellsFolder filesep '' outputShellProjectionName];
            
            writeImage(shellProjection, outputFullPath);
            
            if save3DShells
                outputShellStackName = [inputStackName(1:(end-3)) 'shellStack_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellStackName];
                
                writeImage(shellStack, outputFullPath);
            end;
            
            if saveOrthoMaps
                outputShellOrthoProjectionYZName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionYZ_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionYZName];
                
                writeImage(squeeze(max(shellStack, [], 1)), outputFullPath);
                
                outputShellOrthoProjectionXYaName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXYa_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXYaName];
                
                writeImage(max(shellStack(:, :, 1:round(size(shellStack, 3) / 2)), [], 3), outputFullPath);
                
                outputShellOrthoProjectionXYbName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXYb_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXYbName];
                
                writeImage(max(shellStack(:, :, (round(size(shellStack, 3) / 2) + 1):end), [], 3), outputFullPath);
                
                outputShellOrthoProjectionXZaName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXZa_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXZaName];
                
                writeImage(squeeze(max(shellStack(:, 1:round(size(shellStack, 2) / 2), :), [], 2)), outputFullPath);
                
                outputShellOrthoProjectionXZbName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXZb_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXZbName];
                
                writeImage(squeeze(max(shellStack(:, (round(size(shellStack, 2) / 2) + 1):end, :), [], 2)), outputFullPath);
            end;
            
            %% process child channel(s)
            
            if ~isempty(childChannels)
                for hIndex = 1:numel(childChannels)
                    
                    %% create/load interpolated stack(s) for child channel(s)
                    
                    h = childChannels(hIndex);
                    
                    inputStackName = [inputFileHeader num2str(t, '%.6d') inputFileMiddle num2str(h, '%.2d') inputFileFooter];
                    outputInterpolatedStackName = [inputStackName(1:(end-3)) 'interpolated.klb'];
                    
                    inputFullPath = [inputFolder filesep '' inputHeader num2str(t, '%.6d') inputFooter filesep '' inputStackName];
                    outputFullPath = [isotropicFolder filesep '' outputInterpolatedStackName];
                    
                    if exist(outputFullPath, 'file') == 2
                        disp(['Existing interpolated stack detected for time point ' num2str(t) ', channel ' num2str(h)]);
                        
                        interpolatedStack = readImage(outputFullPath);
                    else
                        stack = readImage(inputFullPath);
                        
                        if showProjections
                            figure; imagesc(max(stack, [], 3));
                            figure; imagesc(squeeze(max(stack, [], 2)));
                        end;
                        
                        if splitFactor == 1
                            xArray = 1:size(stack, 1);
                            yArray = 1:size(stack, 2);
                            zArray = 1:size(stack, 3);
                            zArrayInterpolated = 1:(1/axialScaling):size(stack, 3);
                            
                            xArrayLength = numel(xArray);
                            yArrayLength = numel(yArray);
                            zArrayLength = numel(zArray);
                            zArrayInterpolatedLength = numel(zArrayInterpolated);
                            
                            interpolatedStack = uint16(interp3(...
                                repmat(yArray, xArrayLength, 1, zArrayLength), ...
                                repmat(xArray', 1, yArrayLength, zArrayLength), ...
                                repmat(permute(zArray, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                single(stack), ...
                                repmat(yArray, xArrayLength, 1, zArrayInterpolatedLength), ...
                                repmat(xArray', 1, yArrayLength, zArrayInterpolatedLength), ...
                                repmat(permute(zArrayInterpolated, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                'cubic'));
                        else
                            zArray = 1:size(stack, 3);
                            zArrayInterpolated = 1:(1/axialScaling):size(stack, 3);
                            
                            zArrayLength = numel(zArray);
                            zArrayInterpolatedLength = numel(zArrayInterpolated);
                            
                            interpolatedStack = zeros(size(stack, 1), size(stack, 2), zArrayInterpolatedLength, 'uint16');
                            
                            for x = 1:splitFactor
                                xStart = (x - 1) * ceil(size(stack, 1) / splitFactor) + 1;
                                if x == splitFactor
                                    xStop = size(stack, 1);
                                else
                                    xStop = xStart + ceil(size(stack, 1) / splitFactor) - 1;
                                end;
                                
                                xArray = 1:(xStop - xStart + 1);
                                xArrayLength = numel(xArray);
                                
                                for y = 1:splitFactor
                                    disp(['Processing slab (' num2str(x) ', ' num2str(y) ') for time point ' num2str(t) ', channel ' num2str(h)]);
                                    
                                    yStart = (y - 1) * ceil(size(stack, 2) / splitFactor) + 1;
                                    if y == splitFactor
                                        yStop = size(stack, 2);
                                    else
                                        yStop = yStart + ceil(size(stack, 2) / splitFactor) - 1;
                                    end;
                                    
                                    yArray = 1:(yStop - yStart + 1);
                                    yArrayLength = numel(yArray);
                                    
                                    interpolatedStack(xStart:xStop, yStart:yStop, :) = uint16(interp3(...
                                        repmat(yArray, xArrayLength, 1, zArrayLength), ...
                                        repmat(xArray', 1, yArrayLength, zArrayLength), ...
                                        repmat(permute(zArray, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                        single(stack(xStart:xStop, yStart:yStop, :)), ...
                                        repmat(yArray, xArrayLength, 1, zArrayInterpolatedLength), ...
                                        repmat(xArray', 1, yArrayLength, zArrayInterpolatedLength), ...
                                        repmat(permute(zArrayInterpolated, [3 1 2]), xArrayLength, yArrayLength, 1), ...
                                        'cubic'));
                                end;
                            end;
                        end;
                        
                        writeImage(interpolatedStack, outputFullPath);
                        
                        if showProjections
                            figure; imagesc(max(interpolatedStack, [], 3));
                            figure; imagesc(squeeze(max(interpolatedStack, [], 2)));
                        end;
                    end;
                    
                    %% use shell map(s) to create shell projection(s)/stack(s) for child channel(s)
                    
                    disp(['Creating shell projection for time point ' num2str(t) ', channel ' num2str(h)]);
                    
                    shellProjection = NaN(phiSize, thetaSize);
                    
                    if save3DShells || saveOrthoMaps
                        shellStack = zeros(size(interpolatedStack, 1), size(interpolatedStack, 2), size(interpolatedStack, 3), 'uint16');
                    end;
                    
                    for i = 1:numel(phiRange)
                        for j = 1:numel(thetaRange)
                            if ~isnan(minMap(i, j))
                                rRangeLocal = minMap(i, j):rFineStepSize:maxMap(i, j);
                                
                                switch coordinateSystem
                                    case 0
                                        xArray = round(embryoCenter(1) + rRangeLocal .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                        yArray = round(embryoCenter(2) - rRangeLocal .* cos(thetaRange(j)));
                                        zArray = round(embryoCenter(3) + rRangeLocal .* (sin(thetaRange(j)) * sin(phiRange(i))));
                                    case 1
                                        xArray = round(embryoCenter(1) + rRangeLocal .* (sin(thetaRange(j)) * cos(phiRange(i))));
                                        yArray = round(embryoCenter(2) - rRangeLocal .* (sin(thetaRange(j)) * sin(phiRange(i))));
                                        zArray = round(embryoCenter(3) - rRangeLocal .* cos(thetaRange(j)));
                                end;
                                
                                outOfBoundsFlag = xArray < 1 | xArray > xSize | yArray < 1 | yArray > ySize | zArray < 1 | zArray > zSize;
                                
                                xRegion = xArray(~outOfBoundsFlag);
                                yRegion = yArray(~outOfBoundsFlag);
                                zRegion = zArray(~outOfBoundsFlag);
                                
                                if ~isempty(xRegion)
                                    linearIndices = sub2ind(size(mask), xRegion, yRegion, zRegion);
                                    shellProjection(i, j) = max(interpolatedStack(linearIndices));
                                    
                                    if save3DShells || saveOrthoMaps
                                        shellStack(linearIndices) = interpolatedStack(linearIndices);
                                    end;
                                end;
                            end;
                        end;
                    end;
                    
                    % figure;
                    % imagesc(shellProjection);
                    % xlabel('theta');
                    % ylabel('phi');
                    % set(gca, 'YTick', [1 phiSize]);
                    % set(gca, 'YTickLabel', [phiRange(1) phiRange(phiSize)]);
                    % set(gca, 'XTick', [1 thetaSize]);
                    % set(gca, 'XTickLabel', [thetaRange(1) thetaRange(thetaSize)]);
                    % colorbar;
                    
                    outputShellProjectionName = [inputStackName(1:(end-3)) 'shellProjection_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                    outputFullPath = [shellsFolder filesep '' outputShellProjectionName];
                    
                    writeImage(shellProjection, outputFullPath);
                    
                    if save3DShells
                        outputShellStackName = [inputStackName(1:(end-3)) 'shellStack_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellStackName];
                        
                        writeImage(shellStack, outputFullPath);
                    end;
                    
                    if saveOrthoMaps
                        outputShellOrthoProjectionYZName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionYZ_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionYZName];
                        
                        writeImage(squeeze(max(shellStack, [], 1)), outputFullPath);
                        
                        outputShellOrthoProjectionXYaName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXYa_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXYaName];
                        
                        writeImage(max(shellStack(:, :, 1:round(size(shellStack, 3) / 2)), [], 3), outputFullPath);
                        
                        outputShellOrthoProjectionXYbName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXYb_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXYbName];
                        
                        writeImage(max(shellStack(:, :, (round(size(shellStack, 3) / 2) + 1):end), [], 3), outputFullPath);
                        
                        outputShellOrthoProjectionXZaName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXZa_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXZaName];
                        
                        writeImage(squeeze(max(shellStack(:, 1:round(size(shellStack, 2) / 2), :), [], 2)), outputFullPath);
                        
                        outputShellOrthoProjectionXZbName = [inputStackName(1:(end-3)) 'shellStackOrthoProjectionXZb_' num2str(phiSteps) '_' num2str(thetaSteps) '_' num2str(rShellThickness) medianString smoothString '.klb'];
                        outputFullPath = [shellsFolder filesep '' outputShellOrthoProjectionXZbName];
                        
                        writeImage(squeeze(max(shellStack(:, (round(size(shellStack, 2) / 2) + 1):end, :), [], 2)), outputFullPath);
                    end;
                end;
            end;
        end;
    end;
end;