function calculateDelta(parameterDatabase, t, memoryEstimate)

load(parameterDatabase);

version = 1.02;

configuration = cell(21, 1);

configuration{1}  = version;      configuration{2}  = t;          configuration{3}  = timeClusters; configuration{4}  = dffSampling;     configuration{5} = subOffset;
configuration{6}  = meanFraction; configuration{7}  = scaling;    configuration{8}  = forceZero;    configuration{9}  = referenceOffset;
configuration{10} = medianFlag;   configuration{11} = kernelSize;
configuration{12} = inputString;  configuration{13} = header;     configuration{14} = footer;       configuration{15} = dataType;
configuration{16} = specimen;     configuration{17} = cameras;    configuration{18} = channels;
configuration{19} = inputType;    configuration{20} = outputType;
configuration{21} = [jobMemory(1) memoryEstimate];

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

disp(['Processing time cluster ' num2str(t, '%.4d')]);

if dataType == 0
    referenceStackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timeClusters{t, 1}, '%.6d')];
    referenceStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timeClusters{t, 1}, '%.6d') configurationString '.baseline' inputExtension];
else
    referenceStackPath = [inputString '.registered' filesep header '.TM' num2str(timeClusters{t, 1}, '%.6d') footer];
    referenceStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timeClusters{t, 1}, '%.6d') configurationString '.fusedStack_baseline' inputExtension];
end;

fullReferenceStackPath = [referenceStackPath filesep '' referenceStackName];
referenceStack = double(readImage(fullReferenceStackPath));

if subOffset(1)
    referenceStack = referenceStack - subOffset(2);
end;

for currentTP = timeClusters{t, 2}
    if dataType == 0
        currentStackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(currentTP, '%.6d')];
        currentStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString inputExtension];
        
        outputPath = [inputString '.processed' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(currentTP, '%.6d')];
    else
        currentStackPath = [inputString '.registered' filesep header '.TM' num2str(currentTP, '%.6d') footer];
        currentStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack' inputExtension];
        
        outputPath = [inputString '.processed' filesep header '.TM' num2str(currentTP, '%.6d') footer];
    end;
    
    fullCurrentStackPath = [currentStackPath filesep '' currentStackName];
    currentStack = double(readImage(fullCurrentStackPath));
    
    if subOffset(1)
        currentStack = currentStack - subOffset(2);
    end;
    
    currentStack = uint16(((currentStack - referenceStack) ./ (referenceStack + referenceOffset) + scaling(1)) .* scaling(2));
    if forceZero(1)
        currentStack(referenceStack <= forceZero(2)) = scaling(1) * scaling(2);
    end;
    
    if exist(outputPath, 'dir') ~= 7
        mkdir(outputPath);
    end;
    
    if dataType == 0
        outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString outputExtension];
        outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.xyProjection' outputExtension];
        outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.xzProjection' outputExtension];
        outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.yzProjection' outputExtension];
    else
        outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack' outputExtension];
        outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_xyProjection' outputExtension];
        outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_xzProjection' outputExtension];
        outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_yzProjection' outputExtension];
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
            outputMedianStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.medianFiltered' outputExtension];
            outputMedianXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.medianFiltered_xyProjection' outputExtension];
            outputMedianXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.medianFiltered_xzProjection' outputExtension];
            outputMedianYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.medianFiltered_yzProjection' outputExtension];
        else
            outputMedianStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_medianFiltered' outputExtension];
            outputMedianXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_medianFiltered.xyProjection' outputExtension];
            outputMedianXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_medianFiltered.xzProjection' outputExtension];
            outputMedianYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString '.fusedStack_medianFiltered.yzProjection' outputExtension];
        end;
        fullOutputMedianStackPath = [outputPath filesep '' outputMedianStackName];
        fullOutputMedianXYProjectionPath = [outputPath filesep '' outputMedianXYProjectionName];
        fullOutputMedianXZProjectionPath = [outputPath filesep '' outputMedianXZProjectionName];
        fullOutputMedianYZProjectionPath = [outputPath filesep '' outputMedianYZProjectionName];
        
        writeImage(currentStack, fullOutputMedianStackPath);
        writeImage(max(currentStack, [], 3), fullOutputMedianXYProjectionPath);
        writeImage(squeeze(max(currentStack, [], 2)), fullOutputMedianXZProjectionPath);
        writeImage(squeeze(max(currentStack, [], 1)), fullOutputMedianYZProjectionPath);
    end;
end;

end