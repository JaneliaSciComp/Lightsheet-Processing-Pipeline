inputRoot1 = 'R:\SV1\KM_15-10-26\Mmu_E1_mKate2_20151026_180901.corrected\Results\TimeFused.Corrected';
inputRoot2 = 'R:\SV1\KM_15-10-26\Mmu_E1_mKate2_20151027_121216.corrected\Results\TimeFused.Corrected';

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

outputRoot = 'R:\SV1\KM_15-10-26\Mmu_E1_mKate2_Combined.corrected\Results\TimeFused.Corrected';

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
    
    currentInputFolder = [inputRoot1 '\' inputHeaderA num2str(inputTimepoints1(i), '%.6d') inputFooterA];
    currentOutputFolder = [outputRoot '\' inputHeaderA num2str(outputTimepoints1(i), '%.6d') inputFooterA];
    
    currentInputFilename = [inputHeaderB num2str(inputTimepoints1(i), '%.6d') inputFooterB];
    currentOutputFilename = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') inputFooterB];
    
    currentOutputProjFilenameXY = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterXY];
    currentOutputProjFilenameXZ = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterXZ];
    currentOutputProjFilenameYZ = [inputHeaderB num2str(outputTimepoints1(i), '%.6d') projFooterYZ];
    
    if exist(currentOutputFolder, 'dir') ~= 7
        mkdir(currentOutputFolder);
    end;
    
    if exist([currentOutputFolder '\' currentOutputFilename], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameXY], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameXZ], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameYZ], 'file') ~= 2
        
        stack = readImage([currentInputFolder '\' currentInputFilename]);
        stack = cat(3, zLowSlab, stack, zHighSlab);
        writeImage(stack, [currentOutputFolder '\' currentOutputFilename]);
        
        writeImage(max(stack, [], 3), [currentOutputFolder '\' currentOutputProjFilenameXY]);
        writeImage(squeeze(max(stack, [], 2)), [currentOutputFolder '\' currentOutputProjFilenameXZ]);
        writeImage(squeeze(max(stack, [], 1)), [currentOutputFolder '\' currentOutputProjFilenameYZ]);
    end;
end;

parfor i = 1:numel(inputTimepoints2)
    disp(['Converting time point ' num2str(inputTimepoints2(i)) ' of second time series']);
    
    currentInputFolder = [inputRoot2 '\' inputHeaderA num2str(inputTimepoints2(i), '%.6d') inputFooterA];
    currentOutputFolder = [outputRoot '\' inputHeaderA num2str(outputTimepoints2(i), '%.6d') inputFooterA];
    
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
    
    if exist([currentOutputFolder '\' currentOutputFilename], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameXY], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameXZ], 'file') ~= 2 || ...
            exist([currentOutputFolder '\' currentOutputProjFilenameYZ], 'file') ~= 2
    
        copyfile([currentInputFolder '\' currentInputFilename], [currentOutputFolder '\' currentOutputFilename]);
        copyfile([currentInputFolder '\' currentInputProjFilenameXY], [currentOutputFolder '\' currentOutputProjFilenameXY]);
        copyfile([currentInputFolder '\' currentInputProjFilenameXZ], [currentOutputFolder '\' currentOutputProjFilenameXZ]);
        copyfile([currentInputFolder '\' currentInputProjFilenameYZ], [currentOutputFolder '\' currentOutputProjFilenameYZ]);
    end;
end;

if matlabpool('size') > 0
    matlabpool('close');
end;
disp(' ');