function interpolateStack(parameterDatabase, t, memoryEstimate)

load(parameterDatabase);
timepoint = timepoints(t);

version = 1.01;

configuration = cell(14, 1);

configuration{1}  = version;      configuration{2}  = timepoint;    configuration{3}  = rootFolder;
configuration{4}  = inputHeader1; configuration{5}  = inputHeader2; configuration{6}  = inputHeader3; configuration{7} = inputFooter; configuration{8} = timeStamp;
configuration{9}  = interMode;    configuration{10} = scaling;      configuration{11} = splitting;
configuration{12} = inputType;    configuration{13} = outputType;   configuration{14} = [jobMemory(1) memoryEstimate];        

outputPath = [rootFolder '.Interpolated\' inputHeader1 num2str(timepoint, timeStamp) inputHeader2];
if exist(outputPath, 'dir') ~= 7
    mkdir(outputPath);
end;

save([outputPath '\' inputHeader3 num2str(timepoint, timeStamp) inputFooter '.configuration.mat'], 'configuration');

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

inputName = [rootFolder '\' inputHeader1 num2str(timepoint, timeStamp) inputHeader2 '\' inputHeader3 num2str(timepoint, timeStamp) inputFooter inputExtension];

stack = readImage(inputName);
xSize = size(stack, 1);
ySize = size(stack, 2);
zSize = size(stack, 3);

if interMode(1, 2) == 0
    interpolationKeyword = 'linear';
elseif interMode(1, 2) == 1
    interpolationKeyword = 'spline';
else
    interpolationKeyword = 'cubic';
end;

oldSpacing = 1:zSize;
newSpacing = 1:(1/scaling):zSize;

interpolatedStack = zeros(xSize, ySize, numel(newSpacing), 'uint16');

if interMode(1, 1) == 0
    for x = 1:xSize
        for y = 1:ySize
            interpolatedStack(x, y, :) = uint16(interp1(oldSpacing, double(squeeze(stack(x, y, :))), newSpacing, interpolationKeyword));
        end;
    end;
else
    for i = 1:splitting
        xStart = round((i - 1) * xSize / splitting + 1);
        xStop = round(i * xSize / splitting);
        
        [xOld, yOld, zOld] = meshgrid(1:ySize, 1:(xStop - xStart + 1), oldSpacing);
        [xNew, yNew, zNew] = meshgrid(1:ySize, 1:(xStop - xStart + 1), newSpacing);
        interpolatedStack(xStart:xStop, :, :) = ...
            uint16(interp3(xOld, yOld, zOld, double(stack(xStart:xStop, :, :)), xNew ,yNew, zNew, interpolationKeyword));
    end;
end;

interpolatedXYProjection = max(interpolatedStack, [], 3);
interpolatedXZProjection = squeeze(max(interpolatedStack, [], 2));

outputXYProjectionName = [rootFolder '.Interpolated\' inputHeader1 num2str(timepoint, timeStamp) inputHeader2 '\' inputHeader3 num2str(timepoint, timeStamp) inputFooter '_xyProjection' outputExtension];
outputXZProjectionName = [rootFolder '.Interpolated\' inputHeader1 num2str(timepoint, timeStamp) inputHeader2 '\' inputHeader3 num2str(timepoint, timeStamp) inputFooter '_xzProjection' outputExtension];
outputStackName = [rootFolder '.Interpolated\' inputHeader1 num2str(timepoint, timeStamp) inputHeader2 '\' inputHeader3 num2str(timepoint, timeStamp) inputFooter outputExtension];

writeImage(interpolatedXYProjection, outputXYProjectionName);
writeImage(interpolatedXZProjection, outputXZProjectionName);
writeImage(interpolatedStack, outputStackName);