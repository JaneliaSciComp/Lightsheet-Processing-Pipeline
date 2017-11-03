timepoints       = 0:278;

dataTypes        = {'LargePSF_iter50'; 'SmallPSF_iter20'};
sigmas           = [20 8; 5 2];
kernelFactor     = 3;
splitting        = 5;
writeProjections = 1;

disp(' ');

for t = 1:numel(timepoints)
    timepoint = timepoints(t);
    
    for d = 1:numel(dataTypes)
        sigma1 = sigmas(d, 1);
        sigma2 = sigmas(d, 2);
        
        outputStackName = ['TM' num2str(timepoint, '%.6d') filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00' ...
            '.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_' dataTypes{d} '_lambdaTV000000.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.uint16.klb'];
        
        lastProjectionName = ['ProjectionsFusion' filesep 'Frames' filesep 'Fusion_' dataTypes{d} '.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.TM' num2str(timepoint, '%.3d') '_yz.klb'];
        
        if exist(outputStackName, 'file') ~= 2 || (writeProjections && exist(lastProjectionName, 'file') ~= 2)
            disp(['Missing stacks for data type ' dataTypes{d} ' at time point ' num2str(timepoint)]);
        end;
    end;
end;

disp(' ');