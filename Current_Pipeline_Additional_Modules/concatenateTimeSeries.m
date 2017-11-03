inputRoot1 = 'R:' filesep 'SV1' filesep 'KM_15-10-26' filesep 'Mmu_E1_mKate2_20151026_180901.corrected' filesep 'Results' filesep 'TimeFused.Corrected';
inputRoot2 = 'R:' filesep 'SV1' filesep 'KM_15-10-26' filesep 'Mmu_E1_mKate2_20151027_121216.corrected' filesep 'Results' filesep 'TimeFused.Corrected';

inputTimepoints1 = 0:270;
inputTimepoints2 = 0:378;

inputHeaderA = 'Mmu_E1_mKate2.TM';
inputFooterA = '_timeFused_blending';
inputHeaderB = 'SPM00_TM';
inputFooterB = '_CM00_CM01_CHN00.fusedStack.corrected.klb';

projFooterXY = '_CM00_CM01_CHN00.fusedStack_xyProjection.corrected.klb';
projFooterXZ = '_CM00_CM01_CHN00.fusedStack_xzProjection.corrected.klb';
projFooterYZ = '_CM00_CM01_CHN00.fusedStack_yzProjection.corrected.klb';

zLowExpansion  = 198;
zHighExpansion = 111;

outputRoot = 'R:' filesep 'SV1' filesep 'KM_15-10-26' filesep 'Mmu_E1_mKate2_Combined.corrected' filesep 'Results' filesep 'TimeFused.Corrected';

poolWorkers = 12;

%% main loop

if exist(outputRoot, 'dir') ~= 7
    mkdir(outputRoot);
end;

zLowSlab = zeros(2048, 2048, zLowExpansion, 'uint16');
zHighSlab = zeros(2048, 2048, zHighExpansion, 'uint16');

outputTimepoints1 = inputTimepoints1;
outputTimepoints2 = inputTimepoints2 + inputTimepoints1(end) + 1;

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);
disp(' ');

parfor i = 1:numel(inputTimepoints1)
    disp(['Converting time point ' num2str(inputTimepoints1(i)) ' of first time series']);
    
    currentInputFolder = [inputRoot1 filesep '' inputHeaderA num2str(inputTimepoints1(i), '%.6d') inputFooterA];
    currentOutputFolder = [outputRoot filesep '' inputHeaderA num2str(outputTimepoints1(i), '%.6d') inputFooterA];
    
    currentInputFilename = [inputHeaderB num2str(inputTimepoints1(i), '%.6d') inputFooterB];
    currentOutputFilename = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') inputFooterB];
    
    currentOutputProjFilenameXY = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterXY];
    currentOutputProjFilenameXZ = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterXZ];
    currentOutputProjFilenameYZ = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterYZ];
    
    if exist(currentOutputFolder, 'dir') ~= 7
        mkdir(currentOutputFolder);
    end;
    
    if exist([currentOutputFolder filesep '' currentOutputFilename], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameXY], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameXZ], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameYZ], 'file') ~= 2
        
        stack = readImage([currentInputFolder filesep '' currentInputFilename]);
        stack = cat(3, zLowSlab, stack, zHighSlab);
        writeImage(stack, [currentOutputFolder filesep '' currentOutputFilename]);
        
        writeImage(max(stack, [], 3), [currentOutputFolder filesep '' currentOutputProjFilenameXY]);
        writeImage(squeeze(max(stack, [], 2)), [currentOutputFolder filesep '' currentOutputProjFilenameXZ]);
        writeImage(squeeze(max(stack, [], 1)), [currentOutputFolder filesep '' currentOutputProjFilenameYZ]);
    end;
end;

parfor i = 1:numel(inputTimepoints2)
    disp(['Converting time point ' num2str(inputTimepoints2(i)) ' of second time series']);
    
    currentInputFolder = [inputRoot2 filesep '' inputHeaderA num2str(inputTimepoints2(i), '%.6d') inputFooterA];
    currentOutputFolder = [outputRoot filesep '' inputHeaderA num2str(outputTimepoints2(i), '%.6d') inputFooterA];
    
    currentInputFilename = [inputHeaderB num2str(inputTimepoints2(i), '%.6d') inputFooterB];
    currentOutputFilename = [inputHeaderB num2str(outputTimepoints2(i), '%.6d') inputFooterB];
    
    currentInputProjFilenameXY = [inputHeaderB num2str(inputTimepoints2(i), '%.6d') projFooterXY];
    currentInputProjFilenameXZ = [inputHeaderB num2str(inputTimepoints2(i), '%.6d') projFooterXZ];
    currentInputProjFilenameYZ = [inputHeaderB num2str(inputTimepoints2(i), '%.6d') projFooterYZ];
    
    currentOutputProjFilenameXY = [inputHeaderB num2str(outputTimepoints2(i), '%.6d') projFooterXY];
    currentOutputProjFilenameXZ = [inputHeaderB num2str(outputTimepoints2(i), '%.6d') projFooterXZ];
    currentOutputProjFilenameYZ = [inputHeaderB num2str(outputTimepoints2(i), '%.6d') projFooterYZ];
    
    if exist(currentOutputFolder, 'dir') ~= 7
        mkdir(currentOutputFolder);
    end;
    
    if exist([currentOutputFolder filesep '' currentOutputFilename], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameXY], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameXZ], 'file') ~= 2 || ...
            exist([currentOutputFolder filesep '' currentOutputProjFilenameYZ], 'file') ~= 2
    
        copyfile([currentInputFolder filesep '' currentInputFilename], [currentOutputFolder filesep '' currentOutputFilename]);
        copyfile([currentInputFolder filesep '' currentInputProjFilenameXY], [currentOutputFolder filesep '' currentOutputProjFilenameXY]);
        copyfile([currentInputFolder filesep '' currentInputProjFilenameXZ], [currentOutputFolder filesep '' currentOutputProjFilenameXZ]);
        copyfile([currentInputFolder filesep '' currentInputProjFilenameYZ], [currentOutputFolder filesep '' currentOutputProjFilenameYZ]);
    end;
end;

if matlabpool('size') > 0
    matlabpool('close');
end;
disp(' ');