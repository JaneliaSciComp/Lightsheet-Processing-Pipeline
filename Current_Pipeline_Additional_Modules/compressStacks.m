timepointsToExclude = [];
resamplingFactor = 10;

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'], 'file') ~= 2 || ...
            ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'], 'file') ~= 2 || ...
            ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_LargePSF_iter50_lambdaTV000000.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_SmallPSF_iter20_lambdaTV000000.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        else
            allDoneFlag = ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_LargePSF_iter50_lambdaTV000000.uint16.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_SmallPSF_iter20_lambdaTV000000.uint16.klb'], 'file') == 2;
            if allDoneFlag
                candidates(i) = [];
            end;
        end;
    end;
end;

for i = 1:numel(candidates)
    timepoint = str2num(candidates(i).name(3:end));
    
    if isempty(find(timepointsToExclude == timepoint, 1))
        disp(['converting time point ' num2str(timepoint)]);
        
        for s = 1:10
            switch s
                case 1
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'];
                case 2
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'];
                case 3
                    inputName = [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'];
                case 4
                    inputName = [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.klb'];
                case 5
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'];
                case 6
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'];
                case 7
                    inputName = [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'];
                case 8
                    inputName = [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.klb'];
                case 9
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_LargePSF_iter50_lambdaTV000000.klb'];
                case 10
                    inputName = [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb_dec_LR_multiGPU_MVD_SmallPSF_iter20_lambdaTV000000.klb'];
            end;
            
            stack = readImage(inputName);
            if (max(stack(:)) * resamplingFactor) > (2^16 - 1)
                error(['Resampled intensity range incompatible with uint16 format for stack ' inputName]);
            end;
            
            stack = uint16(stack .* resamplingFactor);
            writeImage(stack, [inputName(1:(end - 3)) 'uint16.klb']);
        end;
    else
        disp(['skipping time point ' num2str(timepoint)]);
    end;
end;