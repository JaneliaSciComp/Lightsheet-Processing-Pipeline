function createMask(parameterDatabase, t, memoryEstimate)

% -----------------------------------------------------------------------------------------------
% | Mask creation                                                                               |
% |                                                                                             |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2017                          |
% | Email: kellerp@janelia.hhmi.org                                                             |
% |                                                                                             |
% | Utilizes optimization modules and functions by Fernando Amat, HHMI/Janelia Research Campus: |
% | readKLBstack.mexw64                                                                         |
% | writeKLBstack.mexw64                                                                        |
% -----------------------------------------------------------------------------------------------

load(parameterDatabase);
timepoint = timepoints(t);

version = 1.00;

configuration = cell(16, 1);

configuration{1}  = version;      configuration{2}  = timepoint;   configuration{3}  = inputDir;    configuration{4}  = outputDir;
configuration{5}  = header;       configuration{6}  = footer;      configuration{7}  = stackLabel;
configuration{8}  = specimen;     configuration{9}  = cameras;     configuration{10} = channels;
configuration{11} = threshold;    configuration{12} = removeDirt;  configuration{13} = connectivity;
configuration{14} = inputType;    configuration{15} = outputType;
configuration{16} = [jobMemory(1) memoryEstimate];        

if length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion)');
end;

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

disp(['processing time point ' num2str(timepoint, '%.6d')]);

fileNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString];
stackName = [inputDir header '.TM' num2str(timepoint, '%.6d') footer filesep '' fileNameHeader '.fusedStack' stackLabel inputExtension];

stack = readImage(stackName);

outputPath = [outputDir header '.TM' num2str(timepoint, '%.6d') footer];
mkdir(outputPath);

save([outputPath filesep '' fileNameHeader '.configuration.mat'], 'configuration');

stackMask = stack > threshold;

if removeDirt
    objectStats = regionprops(bwconncomp(stackMask, 6), 'Area', 'PixelIdxList');
    largestBlobIndex = 1;
    for i = 1:numel(objectStats)
        if objectStats(i).Area > objectStats(largestBlobIndex).Area
            largestBlobIndex = i;
        end;
    end;
    stackMask = false(size(stackMask));
    stackMask(objectStats(largestBlobIndex).PixelIdxList) = 1;
end;

stack(~stackMask) = 0;

xyProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.masked' outputExtension];
xzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.masked' outputExtension];
yzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.masked' outputExtension];

writeImage(max(stack, [], 3), xyProjectedName);
writeImage(squeeze(max(stack, [], 2)), xzProjectedName);
writeImage(squeeze(max(stack, [], 1)), yzProjectedName);

maskStackName = [outputPath filesep '' fileNameHeader '.fusedStack.mask' outputExtension];
writeImage(uint8(stackMask), maskStackName);

xyProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.mask' outputExtension];
xzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.mask' outputExtension];
yzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.mask' outputExtension];

writeImage(uint8(max(stackMask, [], 3)), xyProjectedName);
writeImage(uint8(squeeze(max(stackMask, [], 2))), xzProjectedName);
writeImage(uint8(squeeze(max(stackMask, [], 1))), yzProjectedName);

clear stack stackMask;

end