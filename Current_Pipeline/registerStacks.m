function registerStacks(parameterDatabase, t, memoryEstimate)

load(parameterDatabase);
timepoint = timepoints(t);
medianFlag = medianFlags(t);

version = 1.03;

configuration = cell(16, 1);

configuration{1}  = version;     configuration{2}  = timepoint;   configuration{3}  = references;
configuration{4}  = medianFlag;  configuration{5}  = dffSampling; configuration{6}  = kernelSize;
configuration{7}  = inputString; configuration{8}  = header;      configuration{9}  = footer;     configuration{10} = dataType;
configuration{11} = specimen;    configuration{12} = cameras;     configuration{13} = channels;
configuration{14} = inputType;   configuration{15} = outputType;
configuration{16} = [jobMemory(1) memoryEstimate];        

if length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
elseif length(cameras) == 1 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % single-view data
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion) or 1 channel and 1 camera (single-view data)');
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

disp(['Processing time point ' num2str(timepoint, '%.6d')]);

if dataType == 0
    inputPath = [inputString filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoint, '%.6d')];
    outputPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoint, '%.6d')];
    inputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString inputExtension];
else
    inputPath = [inputString filesep '' header '.TM' num2str(timepoint, '%.6d') footer];
    outputPath = [inputString '.registered' filesep header '.TM' num2str(timepoint, '%.6d') footer];
    inputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack' inputExtension];
end;

if exist(outputPath, 'dir') ~= 7
    mkdir(outputPath);
end;

save([outputPath filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.configuration.mat'], 'configuration');

fullInputStackPath = [inputPath filesep '' inputStackName];
currentStack = readImage(fullInputStackPath);

referencePath = [inputString '.registered' filesep 'Reference'];
referenceStackName = ['SPM' num2str(specimen, '%.2d') configurationString '.referenceStack' outputExtension];
fullReferenceStackPath = [referencePath filesep '' referenceStackName];
referenceStack = readImage(fullReferenceStackPath);

shifts = zeros(size(currentStack, 3), 2);

for z = 1:size(currentStack, 3)
    [dx, dy] = registerImages(referenceStack(:, :, z), currentStack(:, :, z));
    
    if dx >= 0 && dy >= 0
        currentStack(1:(end - dx), 1:(end - dy), z) = currentStack((1 + dx):end, (1 + dy):end, z);
    elseif dx >= 0 && dy < 0
        currentStack(1:(end - dx), (1 - dy):end, z) = currentStack((1 + dx):end, 1:(end + dy), z);
    elseif dx < 0 && dy >= 0
        currentStack((1 - dx):end, 1:(end - dy), z) = currentStack(1:(end + dx), (1 + dy):end, z);
    elseif dx < 0 && dy < 0
        currentStack((1 - dx):end, (1 - dy):end, z) = currentStack(1:(end + dx), 1:(end + dy), z);
    end;
    
    shifts(z, :) = [dx, dy];
end;

if dataType == 0
    outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString outputExtension];
    outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.xyProjection' outputExtension];
    outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.xzProjection' outputExtension];
    outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.yzProjection' outputExtension];
else
    outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack' outputExtension];
    outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack_xyProjection' outputExtension];
    outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack_xzProjection' outputExtension];
    outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack_yzProjection' outputExtension];
end;
fullOutputStackPath = [outputPath filesep '' outputStackName];
fullOutputXYProjectionPath = [outputPath filesep '' outputXYProjectionName];
fullOutputXZProjectionPath = [outputPath filesep '' outputXZProjectionName];
fullOutputYZProjectionPath = [outputPath filesep '' outputYZProjectionName];

writeImage(currentStack, fullOutputStackPath);
writeImage(max(currentStack, [], 3), fullOutputXYProjectionPath);
writeImage(squeeze(max(currentStack, [], 2)), fullOutputXZProjectionPath);
writeImage(squeeze(max(currentStack, [], 1)), fullOutputYZProjectionPath);

if medianFlag
    for z = 1:size(currentStack, 3)
        currentStack(:, :, z) = medfilt2(currentStack(:, :, z), [kernelSize kernelSize]);
    end;
    
    if dataType == 0
        outputMedianStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.medianFiltered' outputExtension];
    else
        outputMedianStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.fusedStack_medianFiltered' outputExtension];
    end;
    fullOutputMedianStackPath = [outputPath filesep '' outputMedianStackName];
    
    writeImage(currentStack, fullOutputMedianStackPath);
end;

save([outputPath filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString '.shifts.mat'], 'shifts');

end

function [dx, dy] = registerImages(im1, im2)

I1 = double(im1);
I1 = I1 - mean(I1(:));
I2 = double(im2);
I2 = I2 - mean(I2(:));

C = ifft2(conj(fft2(I1)) .* (fft2(I2)));
C = fftshift(C);
[x, y] = find(C == max(C(:)), 1);

dx = x - round(size(im1, 1) / 2 - 0.5) - 1;
dy = y - round(size(im1, 2) / 2 - 0.5) - 1;

end