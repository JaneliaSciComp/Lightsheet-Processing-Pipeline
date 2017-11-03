timepoints = 0:278;
forceOverwrite = 0;

roiRootFolder = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ROIs' filesep 'Vectors';

footers = {...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.fusionSigma_20_8.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.fusionSigma_5_2.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_LargePSF_iter50_lambdaTV000000.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'};

projectionFolders = {...
    'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsFusion' filesep 'FramesCropped';...
    'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsFusion' filesep 'FramesCropped';...
    'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsMVD' filesep 'FramesCropped';...
    'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ProjectionsMVD' filesep 'FramesCropped'};

projectionHeaders = {...
    'Fusion_LargePSF_iter50.fusionSigma_20_8.TM';...
    'Fusion_SmallPSF_iter20.fusionSigma_5_2.TM';...
    'MVD_LargePSF_iter50.TM';...
    'MVD_SmallPSF_iter20.TM'};

for i = 1:numel(projectionFolders)
    if exist(projectionFolders{i}, 'dir') ~= 7
        mkdir(projectionFolders{i});
    end;
end;

for t = timepoints
    disp(['cropping time point ' num2str(t)]);
    
    load([roiRootFolder filesep 'TM' num2str(t, '%.6d') '_ROI.mat']);
    
    for i = 1:numel(footers)
        if forceOverwrite || exist([projectionFolders{i} filesep '' projectionHeaders{i} num2str(t, '%.3d') '_xz.klb'], 'file') ~= 2
            stack = readImage(['TM' num2str(t, '%.6d') filesep 'SPM00_TM' num2str(t, '%.6d') footers{i}]);
            stack = stack(croppingVectors(i, 1):croppingVectors(i, 2), croppingVectors(i, 3):croppingVectors(i, 4), croppingVectors(i, 5):croppingVectors(i, 6));
            
            writeImage(stack, ['TM' num2str(t, '%.6d') filesep 'SPM00_TM' num2str(t, '%.6d') footers{i}(1:(end - 3)) 'cropped.klb']);
            writeImage(max(stack, [], 3), [projectionFolders{i} filesep '' projectionHeaders{i} num2str(t, '%.3d') '_xy.klb']);
            writeImage(squeeze(max(stack, [], 2)), [projectionFolders{i} filesep '' projectionHeaders{i} num2str(t, '%.3d') '_xz.klb']);
        end;
    end;
end;