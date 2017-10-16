timepoints = 0:278;

rootFolder = 'V:\SV1\KM_15-08-10\Mmu_E1_mKate2_20150810_160708.corrected.registered';

headers = {...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM01_TM';...
    'SPM01_TM';...
    'SPM01_TM';...
    'SPM01_TM';...
    };

footers = {...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb';...
    '_CM00_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_LargePSF_iter50_lambdaTV000000.uint16.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_SmallPSF_iter20_lambdaTV000000.uint16.klb';...
    };

xArray = zeros(numel(timepoints), numel(headers));
yArray = zeros(numel(timepoints), numel(headers));
zArray = zeros(numel(timepoints), numel(headers));

for t = timepoints
    for i = 1:numel(headers)
        currentHeader = readKLBheader([rootFolder '\TM' num2str(t, '%.6d') '\' headers{i} num2str(t, '%.6d') footers{i}]);
        xArray(find(timepoints == t, 1), i) = currentHeader.xyzct(1);
        yArray(find(timepoints == t, 1), i) = currentHeader.xyzct(2);
        zArray(find(timepoints == t, 1), i) = currentHeader.xyzct(3);
    end;
end;

sizeInconsistencyFlags = zeros(numel(timepoints), 3);
for t = timepoints
    for i = 2:numel(headers)
        if xArray(find(timepoints == t, 1), 1) ~= xArray(find(timepoints == t, 1), i)
            sizeInconsistencyFlags(find(timepoints == t, 1), 1) = 0;
            disp(['inconsistent x size at time point ' num2str(t)]);
            break;
        end;
    end;
    for i = 2:numel(headers)
        if yArray(find(timepoints == t, 1), 1) ~= yArray(find(timepoints == t, 1), i)
            sizeInconsistencyFlags(find(timepoints == t, 1), 2) = 0;
            disp(['inconsistent y size at time point ' num2str(t)]);
            break;
        end;
    end;
    for i = 2:numel(headers)
        if zArray(find(timepoints == t, 1), 1) ~= zArray(find(timepoints == t, 1), i)
            sizeInconsistencyFlags(find(timepoints == t, 1), 3) = 0;
            disp(['inconsistent z size at time point ' num2str(t)]);
            break;
        end;
    end;
end;