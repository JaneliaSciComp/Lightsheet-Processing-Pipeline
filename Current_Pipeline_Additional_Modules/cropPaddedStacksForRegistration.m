timepoints = 0:278;
forceOverwrite = 0;

roiRootFolder = ['V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'ROIs' filesep 'Vectors'];

headers = {...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM01_TM';...
    'SPM01_TM'};

footers = {...
    '_CM00_CHN00.affine.trsf.cropped.padded.klb';...
    '_CM01_CHN00.affine.trsf.cropped.padded.klb';...
    '_CM00_CHN00.affine.trsf.cropped.padded.klb';...
    '_CM01_CHN00.affine.trsf.cropped.padded.klb'};

for t = timepoints
    disp(['cropping time point ' num2str(t)]);
    
    load([roiRootFolder filesep 'TM' num2str(t, '%.6d') '_ROI.mat']);
    
    for i = 1:numel(footers)
        outputName = ['TM' num2str(t, '%.6d') filesep '' headers{i} num2str(t, '%.6d') footers{i}(1:(end - 3)) 'cropped.klb'];
        
        if forceOverwrite || exist(outputName, 'file') ~= 2
            stack = readImage(['TM' num2str(t, '%.6d') filesep '' headers{i} num2str(t, '%.6d') footers{i}]);
            stack = stack(croppingVectors(4, 1):croppingVectors(4, 2), croppingVectors(4, 3):croppingVectors(4, 4), croppingVectors(4, 5):croppingVectors(4, 6));
            writeImage(stack, outputName);
        end;
    end;
end;