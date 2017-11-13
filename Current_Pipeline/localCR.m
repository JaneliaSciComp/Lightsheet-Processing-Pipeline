%% parameters

timepoints   = 0:9000;
dffSampling  = 5;      % time interval in which new dF/F reference stacks are computed
dffRadius    = 20;
percentile   = 30;

inputString  = ['X:' filesep 'SiMView2' filesep '14-01-04' filesep 'Dme_L1_57C10-GCaMP641_0_20140104_114246.corrected' filesep 'Results' filesep 'TimeFused'];
header       = 'Dme_L1_57C10-GCaMP641'; % only required for dataType == 1
footer       = '_timeFused_blending';   % only required for dataType == 1
dataType     = 1;                       % 0: unfused processed data (output of clusterPT), 1: fused processed data (output of clusterMF or clusterTF)
specimen     = 0;
cameras      = [1 0];
channels     = 0;

inputType    = 0;      % 0: input data in KLB format
                       % 1: input data in JP2 format
                       % 2: input data in TIF format
outputType   = 0;      % 0: output data saved in KLB format
                       % 1: output data saved in JP2 format
                       % 2: output data saved in TIF format

poolWorkers  = 0;      % use "0" to enable automated detection of available CPU cores

%% main loop

if poolWorkers == 0
    poolWorkers = feature('numcores');
    disp(' ');
    disp([num2str(poolWorkers) ' CPU cores were detected and will be allocated for parallel processing.']);
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

disp(' ');
disp('Setting up processing environment');

masterClock = tic;

ticks = 1:dffSampling:numel(timepoints);
timeArray = zeros(numel(ticks) + 1, 1);

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);

if dataType == 0
    sampleStackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(1), '%.6d')];
    sampleStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString '.medianFiltered' inputExtension];
else
    sampleStackPath = [inputString '.registered' filesep header '.TM' num2str(timepoints(1), '%.6d') footer];
    sampleStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString '.fusedStack_medianFiltered' inputExtension];
end;
fullSampleStackPath = [sampleStackPath filesep '' sampleStackName];

switch inputType
    case 0
        headerInformation = readKLBheader(fullSampleStackPath);
        stackDimensions = headerInformation.xyzct(1:3);
    case 1
        [stackDimensions, bitDepth] = readJP2header(fullSampleStackPath);
    case 2
        headerInformation = imfinfo(fullSampleStackPath);
        stackDimensions = [headerInformation(1).Height headerInformation(1).Width numel(headerInformation)];
end;

dataPool = NaN(stackDimensions(1), stackDimensions(2), stackDimensions(3), 2 * dffRadius, 'single');

timeArray(1) = toc(masterClock);

disp(['Processing environment set up in ' num2str(timeArray(1), '%.2f') ' seconds']);

for t = 1:numel(ticks)
    disp(' ');
    disp(['Processing baseline data for tick ' num2str(t) ' of ' num2str(numel(ticks))]);
    
    masterClock = tic;
    
    if t == 1
        referenceTimePoints = timepoints(ticks(1:dffRadius));
        
        for r = 1:dffRadius
            currentTime = referenceTimePoints(r);
            disp(['* Loading image data at t = ' num2str(currentTime)]);
            
            if dataType == 0
                stackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(currentTime, '%.6d')];
                stackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTime, '%.6d') configurationString '.medianFiltered' inputExtension];
            else
                stackPath = [inputString '.registered' filesep header '.TM' num2str(currentTime, '%.6d') footer];
                stackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTime, '%.6d') configurationString '.fusedStack_medianFiltered' inputExtension];
            end;
            fullStackPath = [stackPath filesep '' stackName];
            
            dataPool(:, :, :, r) = readImage(fullStackPath);
        end;
        
        activeIndices = 1:dffRadius;
        newSlot = dffRadius + 1;
    else
        newTickIndex = dffRadius + t - 1;
        
        if newTickIndex <= numel(ticks)
            newTime = timepoints(ticks(newTickIndex));
            disp(['* Loading image data at t = ' num2str(newTime)]);
            
            if dataType == 0
                stackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(newTime, '%.6d')];
                stackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(newTime, '%.6d') configurationString '.medianFiltered' inputExtension];
            else
                stackPath = [inputString '.registered' filesep header '.TM' num2str(newTime, '%.6d') footer];
                stackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(newTime, '%.6d') configurationString '.fusedStack_medianFiltered' inputExtension];
            end;
            fullStackPath = [stackPath filesep '' stackName];
            
            dataPool(:, :, :, newSlot) = readImage(fullStackPath);
            
            indexFlag = find(activeIndices == newSlot, 1);
            if isempty(indexFlag)
                activeIndices = cat(2, activeIndices, newSlot);
            end;
        else
            disp(['* Removing image data at slot ' num2str(newSlot)]);
            dataPool(:, :, :, newSlot) = NaN;
            
            indexFlag = find(activeIndices == newSlot, 1);
            if ~isempty(indexFlag)
                activeIndices(indexFlag) = [];
            end;
        end;
        
        newSlot = newSlot + 1;
        if newSlot > (2 * dffRadius)
            newSlot = 1;
        end;
    end;
    
    disp('* Multi-threaded percentile computation');
    referenceStack = zeros(stackDimensions(1), stackDimensions(2), stackDimensions(3), 'uint16');
    parfor z = 1:stackDimensions(3)
        referenceStack(:, :, z) = prctile(dataPool(:, :, z, activeIndices), percentile, 4);
    end;
    
    disp('* Saving reference stack');
    currentTimepoint = timepoints(ticks(t));
    
    if dataType == 0
        outputStackPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(currentTimepoint, '%.6d')];
        outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.baseline' outputExtension];
        outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.baseline_xyProjection' outputExtension];
        outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.baseline_xzProjection' outputExtension];
        outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.baseline_yzProjection' outputExtension];
    else
        outputStackPath = [inputString '.registered' filesep header '.TM' num2str(currentTimepoint, '%.6d') footer];
        outputStackName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.fusedStack_baseline' outputExtension];
        outputXYProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.fusedStack_baseline.xyProjection' outputExtension];
        outputXZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.fusedStack_baseline.xzProjection' outputExtension];
        outputYZProjectionName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTimepoint, '%.6d') configurationString '.fusedStack_baseline.yzProjection' outputExtension];
    end;
    fullOutputStackPath = [outputStackPath filesep '' outputStackName];
    fullOutputXYProjectionPath = [outputStackPath filesep '' outputXYProjectionName];
    fullOutputXZProjectionPath = [outputStackPath filesep '' outputXZProjectionName];
    fullOutputYZProjectionPath = [outputStackPath filesep '' outputYZProjectionName];
    
    writeImage(referenceStack, fullOutputStackPath);
    writeImage(max(referenceStack, [], 3), fullOutputXYProjectionPath);
    writeImage(squeeze(max(referenceStack, [], 2)), fullOutputXZProjectionPath);
    writeImage(squeeze(max(referenceStack, [], 1)), fullOutputYZProjectionPath);
    
    timeArray(t + 1) = toc(masterClock);

    disp(['Tick ' num2str(t) ' of ' num2str(numel(ticks)) ' (time point ' num2str(currentTimepoint) ') processed in ' num2str(timeArray(t + 1), '%.2f') ' seconds']);
end;

disp(' ');

if matlabpool('size') > 0
    matlabpool('close');
end;

disp(' ');

elapsedTime = sum(timeArray);
disp(['Processing completed in ' num2str(elapsedTime / 60, '%.2f') ' minutes']);

disp(' ');