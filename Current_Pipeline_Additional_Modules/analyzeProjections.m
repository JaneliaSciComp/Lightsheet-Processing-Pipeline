projectionFolderFusion = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsFusion' filesep 'Frames';
projectionFolderMVD = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsMVD' filesep 'Frames';

projectionTypesFusion = {...
    'Fusion_LargePSF_iter50.fusionSigma_20_8.TM';...
    'Fusion_SmallPSF_iter20.fusionSigma_5_2.TM'};

projectionTypesMVD = {...
    'MVD_LargePSF_iter50.TM';...
    'MVD_SmallPSF_iter20.TM'};

roiFolder = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered.projections' filesep 'TM';

timepoints        = 0:278;
sigma             = 5;
backgroundPrctile = 10;
threshold         = 0.5;
minSize           = 10^5;
margin            = 40;

if exist('ROIs', 'dir') ~= 7
    mkdir('ROIs');
end;

if exist('ROIs' filesep 'OverlayFrames', 'dir') ~= 7
    mkdir('ROIs' filesep 'OverlayFrames');
end;

if exist('ROIs' filesep 'IndividualFrames', 'dir') ~= 7
    mkdir('ROIs' filesep 'IndividualFrames');
end;

if exist('ROIs' filesep 'Vectors', 'dir') ~= 7
    mkdir('ROIs' filesep 'Vectors');
end;

for t = timepoints    
    if exist(['ROIs' filesep 'Vectors' filesep 'TM' num2str(t, '%.6d') '_ROI.mat'], 'file') ~= 2
        disp(['analyzing projections for time point ' num2str(t)]);
        
        croppingVectors = zeros(numel(projectionTypesFusion) + numel(projectionTypesMVD), 6, 'uint16');
        
        load([roiFolder num2str(t, '%.6d') filesep 'TM' num2str(t, '%.6d') '_ROI.mat']);
        
        %% process xy projections
        
        xyProjections = cell(numel(projectionTypesFusion) + numel(projectionTypesMVD), 1);
        sizeArray = zeros(numel(projectionTypesFusion) + numel(projectionTypesMVD), 2);
        for i = 1:numel(projectionTypesFusion)
            xyProjections{i, 1} = readImage([projectionFolderFusion filesep '' projectionTypesFusion{i} num2str(t, '%.3d') '_xy.klb']);
            xyProjections{i, 1} = xyProjections{i, 1}(1:(croppingVector(2)-croppingVector(1)+1), 1:(croppingVector(4)-croppingVector(3)+1));
            sizeArray(i, :) = [size(xyProjections{i, 1}, 1), size(xyProjections{i, 1}, 2)];
            % figure; imagesc(xyProjections{i, 1});
        end;
        for i = 1:numel(projectionTypesMVD)
            xyProjections{i + numel(projectionTypesFusion), 1} = readImage([projectionFolderMVD filesep '' projectionTypesMVD{i} num2str(t, '%.3d') '_xy.klb']);
            sizeArray(i + numel(projectionTypesFusion), :) = ...
                [size(xyProjections{i + numel(projectionTypesFusion), 1}, 1), size(xyProjections{i + numel(projectionTypesFusion), 1}, 2)];
            % figure; imagesc(xyProjections{i + numel(projectionTypesFusion), 1});
        end;
        minX = min(sizeArray(:, 1));
        minY = min(sizeArray(:, 2));
        equalizingArray = zeros(numel(projectionTypesFusion) + numel(projectionTypesMVD), 4);
        for i = 1:numel(xyProjections)
            equalizingArray(i, 1) = round((sizeArray(i, 1) - minX) / 2) + 1;
            equalizingArray(i, 2) = equalizingArray(i, 1) + minX - 1;
            equalizingArray(i, 3) = round((sizeArray(i, 2) - minY) / 2) + 1;
            equalizingArray(i, 4) = equalizingArray(i, 3) + minY - 1;
            xyProjections{i, 1} = xyProjections{i, 1}(equalizingArray(i, 1):equalizingArray(i, 2), ...
                equalizingArray(i, 3):equalizingArray(i, 4));
            % figure; imagesc(xyProjections{i, 1});
            if i == 1
                sumProjection = single(xyProjections{i, 1});
            else
                sumProjection = sumProjection + single(xyProjections{i, 1});
            end;
        end;
        sumProjection = uint16(sumProjection ./ numel(xyProjections));
        
        sumProjectionBlurred = imgaussianAnisotropy(double(sumProjection), [sigma sigma], [sigma sigma] .* 3);
        meanIntensity = mean(sumProjectionBlurred(:));
        background = prctile(sumProjectionBlurred(:), backgroundPrctile);
        adaptiveTreshold = background + threshold * (meanIntensity - background);
        sumProjectionBlurredThresholded = single(bwareaopen(sumProjectionBlurred > adaptiveTreshold, minSize));
        minX = max(1, find(sum(sumProjectionBlurredThresholded, 2) > 0, 1, 'first') - margin);
        maxX = min(size(sumProjection, 1), find(sum(sumProjectionBlurredThresholded, 2) > 0, 1, 'last') + margin);
        minY = max(1, find(sum(sumProjectionBlurredThresholded, 1) > 0, 1, 'first') - margin);
        maxY = min(size(sumProjection, 2), find(sum(sumProjectionBlurredThresholded, 1) > 0, 1, 'last') + margin);
        xyOverlay = sumProjection(minX:maxX, minY:maxY);
        writeImage(xyOverlay, ['ROIs' filesep 'OverlayFrames' filesep 'TM' num2str(t, '%.6d') '_overlay.xy.cropped.klb']);
        
        for i = 1:numel(xyProjections)
            croppingVectors(i, 1) = equalizingArray(i, 1) + minX - 1;
            croppingVectors(i, 2) = croppingVectors(i, 1) + maxX - minX + 1;
            croppingVectors(i, 3) = equalizingArray(i, 3) + minY - 1;
            croppingVectors(i, 4) = croppingVectors(i, 3) + maxY - minY + 1;
            writeImage(xyProjections{i, 1}(minX:maxX, minY:maxY), ...
                ['ROIs' filesep 'IndividualFrames' filesep 'DataSet' num2str(i) '_TM' num2str(t, '%.6d') '.xy.cropped.klb']);
        end;
        
        minXMaster = minX;
        maxXMaster = maxX;
        
        %% process xz projections (keep only z information)
        
        xzProjections = cell(numel(projectionTypesFusion) + numel(projectionTypesMVD), 1);
        sizeArray = zeros(numel(projectionTypesFusion) + numel(projectionTypesMVD), 2);
        for i = 1:numel(projectionTypesFusion)
            xzProjections{i, 1} = readImage([projectionFolderFusion filesep '' projectionTypesFusion{i} num2str(t, '%.3d') '_xz.klb']);
            xzProjections{i, 1} = xzProjections{i, 1}(1:(croppingVector(2)-croppingVector(1)+1), 1:(croppingVector(6)-croppingVector(5)+1));
            sizeArray(i, :) = [size(xzProjections{i, 1}, 1), size(xzProjections{i, 1}, 2)];
            % figure; imagesc(xzProjections{i, 1});
        end;
        for i = 1:numel(projectionTypesMVD)
            xzProjections{i + numel(projectionTypesFusion), 1} = readImage([projectionFolderMVD filesep '' projectionTypesMVD{i} num2str(t, '%.3d') '_xz.klb']);
            sizeArray(i + numel(projectionTypesFusion), :) = ...
                [size(xzProjections{i + numel(projectionTypesFusion), 1}, 1), size(xzProjections{i + numel(projectionTypesFusion), 1}, 2)];
            % figure; imagesc(xzProjections{i + numel(projectionTypesFusion), 1});
        end;
        minX = min(sizeArray(:, 1));
        minZ = min(sizeArray(:, 2));
        equalizingArray = zeros(numel(projectionTypesFusion) + numel(projectionTypesMVD), 4);
        for i = 1:numel(xzProjections)
            equalizingArray(i, 1) = round((sizeArray(i, 1) - minX) / 2) + 1;
            equalizingArray(i, 2) = equalizingArray(i, 1) + minX - 1;
            equalizingArray(i, 3) = round((sizeArray(i, 2) - minZ) / 2) + 1;
            equalizingArray(i, 4) = equalizingArray(i, 3) + minZ - 1;
            xzProjections{i, 1} = xzProjections{i, 1}(equalizingArray(i, 1):equalizingArray(i, 2), ...
                equalizingArray(i, 3):equalizingArray(i, 4));
            % figure; imagesc(xzProjections{i, 1});
            if i == 1
                sumProjection = single(xzProjections{i, 1});
            else
                sumProjection = sumProjection + single(xzProjections{i, 1});
            end;
        end;
        sumProjection = uint16(sumProjection ./ numel(xzProjections));
        
        sumProjectionBlurred = imgaussianAnisotropy(double(sumProjection), [sigma sigma], [sigma sigma] .* 3);
        meanIntensity = mean(sumProjectionBlurred(:));
        background = prctile(sumProjectionBlurred(:), backgroundPrctile);
        adaptiveTreshold = background + threshold * (meanIntensity - background);
        sumProjectionBlurredThresholded = single(bwareaopen(sumProjectionBlurred > adaptiveTreshold, minSize));
        minX = max(1, find(sum(sumProjectionBlurredThresholded, 2) > 0, 1, 'first') - margin);
        maxX = min(size(sumProjection, 1), find(sum(sumProjectionBlurredThresholded, 2) > 0, 1, 'last') + margin);
        minZ = max(1, find(sum(sumProjectionBlurredThresholded, 1) > 0, 1, 'first') - margin);
        maxZ = min(size(sumProjection, 2), find(sum(sumProjectionBlurredThresholded, 1) > 0, 1, 'last') + margin);
        xzOverlay = sumProjection(minX:maxX, minZ:maxZ);
        writeImage(xzOverlay, ['ROIs' filesep 'OverlayFrames' filesep 'TM' num2str(t, '%.6d') '_overlay.xz.cropped.klb']);
        
        for i = 1:numel(xzProjections)
            croppingVectors(i, 5) = equalizingArray(i, 3) + minZ - 1;
            croppingVectors(i, 6) = croppingVectors(i, 5) + maxZ - minZ + 1;
            writeImage(xzProjections{i, 1}(minXMaster:maxXMaster, minZ:maxZ), ...
                ['ROIs' filesep 'IndividualFrames' filesep 'DataSet' num2str(i) '_TM' num2str(t, '%.6d') '.xz.cropped.klb']);
        end;
        
        save(['ROIs' filesep 'Vectors' filesep 'TM' num2str(t, '%.6d') '_ROI.mat'], 'croppingVectors');
    end;
end;