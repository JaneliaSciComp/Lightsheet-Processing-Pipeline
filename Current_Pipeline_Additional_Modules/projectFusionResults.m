timepoints = 0:278;

frames = cell(numel(timepoints), 6);
sizeArray = zeros(numel(timepoints), 12);
intensity = zeros(numel(timepoints), 6, 2);
intensity(:, 1, 1) = timepoints';
intensity(:, 1, 2) = timepoints';

for i = 1:numel(timepoints)
    timepoint = timepoints(i);
    
    disp(['reading time point ' num2str(timepoint)]);
    
    frames{i, 1} = readImage(['ProjectionsFusion\Frames\Fusion_SmallPSF_iter20.fusionSigma_5_2.TM' num2str(timepoint, '%.3d') '_xy.klb']);
    frames{i, 2} = readImage(['ProjectionsFusion\Frames\Fusion_SmallPSF_iter20.fusionSigma_5_2.TM' num2str(timepoint, '%.3d') '_xz.klb']);
    frames{i, 3} = readImage(['ProjectionsFusion\Frames\Fusion_SmallPSF_iter20.fusionSigma_5_2.TM' num2str(timepoint, '%.3d') '_yz.klb']);
    
    sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
    sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
    sizeArray(i, 5:6) = [size(frames{i, 3}, 1) size(frames{i, 3}, 2)];
    
    intensity(i, 2:6, 1) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
    
    frames{i, 4} = readImage(['ProjectionsFusion\Frames\Fusion_LargePSF_iter50.fusionSigma_20_8.TM' num2str(timepoint, '%.3d') '_xy.klb']);
    frames{i, 5} = readImage(['ProjectionsFusion\Frames\Fusion_LargePSF_iter50.fusionSigma_20_8.TM' num2str(timepoint, '%.3d') '_xz.klb']);
    frames{i, 6} = readImage(['ProjectionsFusion\Frames\Fusion_LargePSF_iter50.fusionSigma_20_8.TM' num2str(timepoint, '%.3d') '_yz.klb']);
    
    sizeArray(i, 7:8) = [size(frames{i, 4}, 1) size(frames{i, 4}, 2)];
    sizeArray(i, 9:10) = [size(frames{i, 5}, 1) size(frames{i, 5}, 2)];
    sizeArray(i, 11:12) = [size(frames{i, 6}, 1) size(frames{i, 6}, 2)];
    
    intensity(i, 2:6, 2) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
end;

save('ProjectionsFusion\IntensityStatistics.mat', 'intensity');

disp('writing projection overlay stacks');

for c = 1:6
    aMax = max(sizeArray(:, (c-1)*2+1));
    bMax = max(sizeArray(:, (c-1)*2+2));
    
    projections = zeros(aMax, bMax, numel(timepoints), 'uint16');
    
    for i = 1:numel(timepoints)
        projections(...
            (round((aMax - sizeArray(i, (c-1)*2+1)) / 2) + 1):(round((aMax - sizeArray(i, (c-1)*2+1)) / 2) + sizeArray(i, (c-1)*2+1)), ...
            (round((bMax - sizeArray(i, (c-1)*2+2)) / 2) + 1):(round((bMax - sizeArray(i, (c-1)*2+2)) / 2) + sizeArray(i, (c-1)*2+2)), ...
            i) = frames{i, c};
    end;
    
    switch c
        case 1
            outputName = 'ProjectionsFusion\Fusion_SmallPSF_iter20.overlay_xy';
        case 2
            outputName = 'ProjectionsFusion\Fusion_SmallPSF_iter20.overlay_xz';
        case 3
            outputName = 'ProjectionsFusion\Fusion_SmallPSF_iter20.overlay_yz';
        case 4
            outputName = 'ProjectionsFusion\Fusion_LargePSF_iter50.overlay_xy';
        case 5
            outputName = 'ProjectionsFusion\Fusion_LargePSF_iter50.overlay_xz';
        case 6
            outputName = 'ProjectionsFusion\Fusion_LargePSF_iter50.overlay_yz';
    end;
    
    writeImage(projections, [outputName '.klb']);
end;