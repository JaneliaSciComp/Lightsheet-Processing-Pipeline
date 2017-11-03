roiRootFolder = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered.projections';
timepointsToExclude = [];

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        else
            allDoneFlag = ...
                exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'], 'file') == 2 && ...
                exist([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'], 'file') == 2 && ...
                exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'], 'file') == 2 && ...
                exist([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'], 'file') == 2;
            if allDoneFlag
                candidates(i) = [];
            end;
        end;
    end;
end;

for i = 1:numel(candidates)
    timepoint = str2num(candidates(i).name(3:end));
    
    if isempty(find(timepointsToExclude == timepoint, 1))
        disp(['cropping time point ' num2str(timepoint)]);
        
        load([roiRootFolder filesep 'TM' num2str(timepoint, '%.6d') filesep 'TM' num2str(timepoint, '%.6d') '_ROI.mat'], 'croppingVector');
        
        stack = readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.klb']);
        stack = stack(croppingVector(1):croppingVector(2), croppingVector(3):croppingVector(4), croppingVector(5):croppingVector(6));
        writeImage(stack, [candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb']);
        
        stack = readImage([candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.klb']);
        stack = stack(croppingVector(1):croppingVector(2), croppingVector(3):croppingVector(4), croppingVector(5):croppingVector(6));
        writeImage(stack, [candidates(i).name filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb']);
        
        stack = readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.klb']);
        stack = stack(croppingVector(1):croppingVector(2), croppingVector(3):croppingVector(4), croppingVector(5):croppingVector(6));
        writeImage(stack, [candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb']);
        
        stack = readImage([candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.klb']);
        stack = stack(croppingVector(1):croppingVector(2), croppingVector(3):croppingVector(4), croppingVector(5):croppingVector(6));
        writeImage(stack, [candidates(i).name filesep 'SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb']);
    else
        disp(['skipping time point ' num2str(timepoint)]);
    end;
end;