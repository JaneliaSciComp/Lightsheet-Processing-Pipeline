timepoints = 0:278;

if exist('ProjectionsMVD', 'dir') ~= 7
    mkdir('ProjectionsMVD');
end;

if exist('ProjectionsMVD' filesep 'Frames', 'dir') ~= 7
    mkdir('ProjectionsMVD' filesep 'Frames');
end;

frames = cell(numel(timepoints), 6);
sizeArray = zeros(numel(timepoints), 12);
intensity = zeros(numel(timepoints), 6, 2);
intensity(:, 1, 1) = timepoints';
intensity(:, 1, 2) = timepoints';

for i = 1:numel(timepoints)
    timepoint = timepoints(i);
    
    disp(['reading time point ' num2str(timepoint)]);
    
    allDone = ...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xy.klb'], 'file') == 2 &&...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xz.klb'], 'file') == 2 &&...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_yz.klb'], 'file') == 2;
    
    if allDone
        frames{i, 1} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xy.klb']);
        frames{i, 2} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xz.klb']);
        frames{i, 3} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_yz.klb']);
        
        sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
        sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
        sizeArray(i, 5:6) = [size(frames{i, 3}, 1) size(frames{i, 3}, 2)];
    else
        stack = readImage(['TM' num2str(timepoint, '%.6d') filesep 'SPM00_TM' num2str(timepoint, '%.6d') ...
            '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_SmallPSF_iter20_lambdaTV000000.uint16.klb']);
        
        frames{i, 1} = max(stack, [], 3);
        frames{i, 2} = squeeze(max(stack, [], 2));
        frames{i, 3} = squeeze(max(stack, [], 1));
        
        sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
        sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
        sizeArray(i, 5:6) = [size(frames{i, 3}, 1) size(frames{i, 3}, 2)];
        
        writeImage(frames{i, 1}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xy.klb']);
        writeImage(frames{i, 2}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_xz.klb']);
        writeImage(frames{i, 3}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_yz.klb']);
    end;
    
    intensity(i, 2:6, 1) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
    
    allDone = ...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xy.klb'], 'file') == 2 &&...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xz.klb'], 'file') == 2 &&...
        exist(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_yz.klb'], 'file') == 2;
    
    if allDone
        frames{i, 4} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xy.klb']);
        frames{i, 5} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xz.klb']);
        frames{i, 6} = readImage(['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_yz.klb']);
        
        sizeArray(i, 7:8) = [size(frames{i, 4}, 1) size(frames{i, 4}, 2)];
        sizeArray(i, 9:10) = [size(frames{i, 5}, 1) size(frames{i, 5}, 2)];
        sizeArray(i, 11:12) = [size(frames{i, 6}, 1) size(frames{i, 6}, 2)];
    else
        stack = readImage(['TM' num2str(timepoint, '%.6d') filesep 'SPM00_TM' num2str(timepoint, '%.6d') ...
            '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_LargePSF_iter50_lambdaTV000000.uint16.klb']);
        
        frames{i, 4} = max(stack, [], 3);
        frames{i, 5} = squeeze(max(stack, [], 2));
        frames{i, 6} = squeeze(max(stack, [], 1));
        
        sizeArray(i, 7:8) = [size(frames{i, 4}, 1) size(frames{i, 4}, 2)];
        sizeArray(i, 9:10) = [size(frames{i, 5}, 1) size(frames{i, 5}, 2)];
        sizeArray(i, 11:12) = [size(frames{i, 6}, 1) size(frames{i, 6}, 2)];
        
        writeImage(frames{i, 4}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xy.klb']);
        writeImage(frames{i, 5}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_xz.klb']);
        writeImage(frames{i, 6}, ['ProjectionsMVD' filesep 'Frames' filesep 'MVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_yz.klb']);
    end;
    
    intensity(i, 2:6, 2) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
end;

save('ProjectionsMVD' filesep 'IntensityStatistics.mat', 'intensity');

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
            outputName = 'ProjectionsMVD' filesep 'MVD_SmallPSF_iter20.overlay_xy';
        case 2
            outputName = 'ProjectionsMVD' filesep 'MVD_SmallPSF_iter20.overlay_xz';
        case 3
            outputName = 'ProjectionsMVD' filesep 'MVD_SmallPSF_iter20.overlay_yz';
        case 4
            outputName = 'ProjectionsMVD' filesep 'MVD_LargePSF_iter50.overlay_xy';
        case 5
            outputName = 'ProjectionsMVD' filesep 'MVD_LargePSF_iter50.overlay_xz';
        case 6
            outputName = 'ProjectionsMVD' filesep 'MVD_LargePSF_iter50.overlay_yz';
    end;
    
    writeImage(projections, [outputName '.klb']);
end;