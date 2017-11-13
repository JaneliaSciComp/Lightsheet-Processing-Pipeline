timepoints = 0:212; % 0:278;

if exist('ProjectionsSVD', 'dir') ~= 7
    mkdir('ProjectionsSVD');
end;

if exist(['ProjectionsSVD' filesep 'FramesTIF'], 'dir') ~= 7
    mkdir(['ProjectionsSVD' filesep 'FramesTIF']);
end;

if exist(['ProjectionsSVD' filesep 'FramesKLB'], 'dir') ~= 7
    mkdir(['ProjectionsSVD' filesep 'FramesKLB']);
end;

for s = 1:4
    switch s
        case 1
            specimenString = 'SPM00';
            cameraString = 'CM00';
        case 2
            specimenString = 'SPM00';
            cameraString = 'CM01';
        case 3
            specimenString = 'SPM01';
            cameraString = 'CM00';
        case 4
            specimenString = 'SPM01';
            cameraString = 'CM01';
    end;            
    
    disp(' ');
    disp(['processing ' specimenString ' ' cameraString]);
    
    frames = cell(numel(timepoints), 6);
    sizeArray = zeros(numel(timepoints), 12);
    intensity = zeros(numel(timepoints), 6, 2);
    intensity(:, 1, 1) = timepoints';
    intensity(:, 1, 2) = timepoints';
    
    for i = 1:numel(timepoints)
        timepoint = timepoints(i);
        
        disp(['reading time point ' num2str(timepoint)]);
        
        allDone = ...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb'], 'file') == 2 &&...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb'], 'file') == 2 &&...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb'], 'file') == 2;
        
        if allDone
            frames{i, 1} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb']);
            frames{i, 2} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb']);
            frames{i, 3} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb']);
            
            sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
            sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
            sizeArray(i, 5:6) = [size(frames{i, 3}, 1) size(frames{i, 3}, 2)];
        else        
            stack = readImage(['TM' num2str(timepoint, '%.6d') filesep '' specimenString '_TM' num2str(timepoint, '%.6d') ...
                '_' cameraString '_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb']);
            
            frames{i, 1} = max(stack, [], 3);
            frames{i, 2} = squeeze(max(stack, [], 2));
            frames{i, 3} = squeeze(max(stack, [], 1));
            
            sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
            sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
            sizeArray(i, 5:6) = [size(frames{i, 3}, 1) size(frames{i, 3}, 2)];
            
            writeImage(frames{i, 1}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.tif']);
            writeImage(frames{i, 2}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.tif']);
            writeImage(frames{i, 3}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.tif']);
            
            writeImage(frames{i, 1}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb']);
            writeImage(frames{i, 2}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb']);
            writeImage(frames{i, 3}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_SmallPSF_iter20.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb']);
        end;
        
        intensity(i, 2:6, 1) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
        
        allDone = ...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb'], 'file') == 2 &&...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb'], 'file') == 2 &&...
            exist(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb'], 'file') == 2;
        
        if allDone
            frames{i, 4} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb']);
            frames{i, 5} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb']);
            frames{i, 6} = readImage(['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb']);
            
            sizeArray(i, 7:8) = [size(frames{i, 4}, 1) size(frames{i, 4}, 2)];
            sizeArray(i, 9:10) = [size(frames{i, 5}, 1) size(frames{i, 5}, 2)];
            sizeArray(i, 11:12) = [size(frames{i, 6}, 1) size(frames{i, 6}, 2)];
        else        
            stack = readImage(['TM' num2str(timepoint, '%.6d') filesep '' specimenString '_TM' num2str(timepoint, '%.6d') ...
                '_' cameraString '_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb']);
            
            frames{i, 4} = max(stack, [], 3);
            frames{i, 5} = squeeze(max(stack, [], 2));
            frames{i, 6} = squeeze(max(stack, [], 1));
            
            sizeArray(i, 7:8) = [size(frames{i, 4}, 1) size(frames{i, 4}, 2)];
            sizeArray(i, 9:10) = [size(frames{i, 5}, 1) size(frames{i, 5}, 2)];
            sizeArray(i, 11:12) = [size(frames{i, 6}, 1) size(frames{i, 6}, 2)];
            
            writeImage(frames{i, 4}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.tif']);
            writeImage(frames{i, 5}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.tif']);
            writeImage(frames{i, 6}, ['ProjectionsSVD' filesep 'FramesTIF' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.tif']);
            
            writeImage(frames{i, 4}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xy.klb']);
            writeImage(frames{i, 5}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_xz.klb']);
            writeImage(frames{i, 6}, ['ProjectionsSVD' filesep 'FramesKLB' filesep 'SVD_LargePSF_iter50.TM' num2str(timepoint, '%.3d') '_' specimenString '_' cameraString '_yz.klb']);
        end;
        
        intensity(i, 2:6, 2) = [min(frames{i, 1}(:)), prctile(frames{i, 1}(:), 10), median(frames{i, 1}(:)), prctile(frames{i, 1}(:), 99.9), max(frames{i, 1}(:))];
    end;
    
    save(['ProjectionsSVD' filesep 'IntensityStatistics.' specimenString '_' cameraString '.mat'], 'intensity');
    
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
                outputName = ['ProjectionsSVD' filesep 'SVD_SmallPSF_iter20.overlay_' specimenString '_' cameraString '_xy'];
            case 2
                outputName = ['ProjectionsSVD' filesep 'SVD_SmallPSF_iter20.overlay_' specimenString '_' cameraString '_xz'];
            case 3
                outputName = ['ProjectionsSVD' filesep 'SVD_SmallPSF_iter20.overlay_' specimenString '_' cameraString '_yz'];
            case 4
                outputName = ['ProjectionsSVD' filesep 'SVD_LargePSF_iter50.overlay_' specimenString '_' cameraString '_xy'];
            case 5
                outputName = ['ProjectionsSVD' filesep 'SVD_LargePSF_iter50.overlay_' specimenString '_' cameraString '_xz'];
            case 6
                outputName = ['ProjectionsSVD' filesep 'SVD_LargePSF_iter50.overlay_' specimenString '_' cameraString '_yz'];
        end;
        
        writeImage(projections, [outputName '.klb']);
        writeImage(projections, [outputName '.tif']);
    end;
end;

disp(' ');