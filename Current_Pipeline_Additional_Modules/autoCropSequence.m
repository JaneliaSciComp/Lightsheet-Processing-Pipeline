sigma             = 5;
backgroundPrctile = 10;
threshold         = 0.5;
minSize           = 10^5;
margin            = 60;

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..') || strcmp(candidates(i).name, 'Overlays')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        else
            doneFlag = ...
                exist([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_ROI.mat'], 'file') == 2;
            if doneFlag
                candidates(i) = [];
            end;
        end;
    end;
end;

for i = 1:numel(candidates)
    timepoint = str2num(candidates(i).name(3:end));
    
    disp(['processing time point ' num2str(timepoint)]);
    
    croppingVector = zeros(1, 6, 'uint16');
    
    frame = single(readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb']));
    xyOverlay = uint16(frame ./ 4);
    writeImage(xyOverlay, [candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.xy.klb']);
    
    xyOverlayBlurred = imgaussianAnisotropy(double(xyOverlay), [sigma sigma], [sigma sigma] .* 3);
    meanIntensity = mean(xyOverlayBlurred(:));
    background = prctile(xyOverlayBlurred(:), backgroundPrctile);
    adaptiveTreshold = background + threshold * (meanIntensity - background);
    xyOverlayBlurredThresholded = single(bwareaopen(xyOverlayBlurred > adaptiveTreshold, minSize));
    minX = max(1, find(sum(xyOverlayBlurredThresholded, 2) > 0, 1, 'first') - margin);
    maxX = min(size(frame, 1), find(sum(xyOverlayBlurredThresholded, 2) > 0, 1, 'last') + margin);
    minY = max(1, find(sum(xyOverlayBlurredThresholded, 1) > 0, 1, 'first') - margin);
    maxY = min(size(frame, 2), find(sum(xyOverlayBlurredThresholded, 1) > 0, 1, 'last') + margin);
    xyOverlay = xyOverlay(minX:maxX, minY:maxY);
    writeImage(xyOverlay, [candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.xy.cropped.klb']);
    
    croppingVector(1:4) = [minX maxX minY maxY];
    
    frame = single(readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb']));
    yzOverlay = uint16(frame ./ 4);
    writeImage(yzOverlay, [candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.yz.klb']);
    
    yzOverlayBlurred = imgaussianAnisotropy(double(yzOverlay), [sigma sigma], [sigma sigma] .* 3);
    meanIntensity = mean(yzOverlayBlurred(:));
    background = prctile(yzOverlayBlurred(:), backgroundPrctile);
    adaptiveTreshold = background + threshold * (meanIntensity - background);
    yzOverlayBlurredThresholded = single(bwareaopen(yzOverlayBlurred > adaptiveTreshold, minSize));
    minY = max(1, find(sum(yzOverlayBlurredThresholded, 2) > 0, 1, 'first') - margin);
    maxY = min(size(frame, 1), find(sum(yzOverlayBlurredThresholded, 2) > 0, 1, 'last') + margin);
    minZ = max(1, find(sum(yzOverlayBlurredThresholded, 1) > 0, 1, 'first') - margin);
    maxZ = min(size(frame, 2), find(sum(yzOverlayBlurredThresholded, 1) > 0, 1, 'last') + margin);
    yzOverlay = yzOverlay(minY:maxY, minZ:maxZ);
    writeImage(yzOverlay, [candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.yz.cropped.klb']);
    
    croppingVector(5:6) = [minZ maxZ];
    
    save([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_ROI.mat'], 'croppingVector');
end;