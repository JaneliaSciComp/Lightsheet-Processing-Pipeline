%% parameters

timepoints   = 0:1458;
dffSampling  = 5;        % time interval in which new dF/F reference stacks are computed
subOffset    = [0 100];  % if first value is 1, second value is subtracted from all stacks prior to further processing
                         % note: only required for unfused image stacks that have not been corrected for camera dark current (typically ~100 grey levels for Orca Flash 4.0)
meanFraction = 0.1;      % fraction of mean reference stack intensity that is added to the denominator in the dF/F calculation
scaling      = [3 5000]; % offset (slot 1) and multiplier (slot 2) used to adjust the dynamic range of the dF/F calculation
forceZero    = [1 20];   % if first value is 1, overwrites signal stack voxels with zeros if the reference stack has intensity
                         % at the respective location lower than or equal to the value provided in the second slot (safety feature for segmented data)
medianFlag   = 0;        % flag for generating and saving median-filtered dF/F stacks and their projections
kernelSize   = 3;        % kernel size of median filter, only required if medianFlag == 1

inputString  = 'X:' filesep 'SiMView2' filesep '14-01-04' filesep 'Dme_L1_57C10-GCaMP641_0_20140104_114246.corrected' filesep 'Results' filesep 'TimeFused';
header       = 'Dme_L1_57C10-GCaMP641'; % only required for dataType == 1
footer       = '_timeFused_blending';   % only required for dataType == 1
dataType     = 1;                       % 0: unfused processed data (output of clusterPT), 1: fused processed data (output of clusterMF or clusterTF)
specimen     = 0;
cameras      = [1 0];
channels     = 0;

inputType    = 1;      % 0: input data in KLB format
                       % 1: input data in JP2 format
                       % 2: input data in TIF format
outputType   = 1;      % 0: output data saved in KLB format
                       % 1: output data saved in JP2 format
                       % 2: output data saved in TIF format

localRun     = [0 0];  % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                       % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                       %         note: use "0" to enable automated detection of available CPU cores

jobMemory    = [1 0];  % slot 1: flag for automated memory management (0: disable, 1: enable)
                       % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                       %         note 1: slot 2 is only evaluated if automated memory management is disabled
                       %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)

%% job submission

if localRun(1) == 1 && localRun(2) == 0
    localRun(2) = feature('numcores');
    disp(' ');
    disp([num2str(localRun(2)) ' CPU cores were detected and will be allocated for parallel processing.']);
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

allTimepoints = timepoints;
ticks = 1:dffSampling:numel(timepoints);
existingResultsDetected = 0;

for currentTP = length(timepoints):-1:1
    if dataType == 0
        outputPath = [inputString '.processed' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(currentTP), '%.6d')];
        if medianFlag == 0
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString '.yzProjection' outputExtension];
        else
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString '.medianFiltered_yzProjection' outputExtension];
        end;
    else
        outputPath = [inputString '.processed' filesep header '.TM' num2str(timepoints(currentTP), '%.6d') footer];
        if medianFlag == 0
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString '.fusedStack_yzProjection' outputExtension];
        else
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString '.fusedStack_medianFiltered.yzProjection' outputExtension];
        end;
    end;
    fullFilePath = [outputPath filesep '' fileName];
    
    if exist(fullFilePath, 'file') == 2
        timepoints(currentTP) = [];
        existingResultsDetected = 1;
    end;
end;

if ~isempty(timepoints)
    nAllTicks = numel(ticks);
    timeClusters = cell(nAllTicks, 2);
    
    for n = nAllTicks:-1:1
        startTime = max(allTimepoints(1), allTimepoints(ticks(n)) - floor(dffSampling / 2));
        if n == nAllTicks
            stopTime = allTimepoints(end);
        else
            stopTime = min(allTimepoints(end), allTimepoints(ticks(n)) + ceil(dffSampling / 2) - 1); 
        end;
        
        currentCluster = startTime:stopTime;
        for c = numel(currentCluster):-1:1
            if isempty(find(timepoints == currentCluster(c), 1))
                currentCluster(c) = [];
            end;
        end;
        
        if ~isempty(currentCluster)
            timeClusters{n, 1} = allTimepoints(ticks(n));
            timeClusters{n, 2} = currentCluster;
        else
            timeClusters(n, :) = [];
        end;
    end;
    
    nTicks = size(timeClusters, 1);
    
    if ~existingResultsDetected
        referenceTime = allTimepoints(ticks(round(numel(ticks) / 2)));
        
        disp(' ');
        disp(['Determining denominator offset for dF/F calculation using reference stack at time point ' num2str(referenceTime) '.']);
        
        outputPath = [inputString '.processed'];
        if exist(outputPath, 'dir') == 7
            warning('Default output folder already exists. Abort processing if existing data should not be overwritten.');
        else
            mkdir(outputPath);
        end;
        
        if dataType == 0
            inputPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(referenceTime, '%.6d')];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(referenceTime, '%.6d') configurationString '.baseline' inputExtension];
        else
            inputPath = [inputString '.registered' filesep header '.TM' num2str(referenceTime, '%.6d') footer];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(referenceTime, '%.6d') configurationString '.fusedStack_baseline' inputExtension];
        end;
        fullFilePath = [inputPath filesep '' fileName];
        
        currentStack = readImage(fullFilePath);
        if subOffset(1)
            currentStack = currentStack - subOffset(2);
        end;
        referenceOffset = ceil(meanFraction * mean(currentStack(currentStack > 0)));
        
        referencePath = [outputPath filesep 'Reference'];
        if exist(referencePath, 'dir') ~= 7
            mkdir(referencePath);
        end;
        
        save([referencePath filesep 'SPM' num2str(specimen, '%.2d') configurationString '.referenceOffset.mat'], 'referenceOffset');
    else
        disp(' ');
        disp('Existing processing results detected. Skipping calculation of denominator offset.');
        
        referencePath = [inputString '.processed' filesep 'Reference'];
        load([referencePath filesep 'SPM' num2str(specimen, '%.2d') configurationString '.referenceOffset.mat']);
    end;
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        if dataType == 0
            inputPath = [inputString '.registered' filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(1), '%.6d')];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString inputExtension];
        else
            inputPath = [inputString '.registered' filesep header '.TM' num2str(timepoints(1), '%.6d') footer];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString '.fusedStack' inputExtension];
        end;
        fullFilePath = [inputPath filesep '' fileName];
        
        try
            switch inputType
                case 0
                    headerInformation = readKLBheader(fullFilePath);
                    stackDimensions = headerInformation.xyzct(1:3);
                case 1
                    [stackDimensions, bitDepth] = readJP2header(fullFilePath);
                case 2
                    headerInformation = imfinfo(fullFilePath);
                    stackDimensions = [headerInformation(1).Height headerInformation(1).Width numel(headerInformation)];
            end;
            unitX = 2 * stackDimensions(1) * stackDimensions(2) * stackDimensions(3) / (1024 ^ 3);
        catch errorMessage
            switch inputType
                case 0
                    error('Failed to open KLB file.');
                case 1
                    error('Failed to open JP2 file.');
                case 2
                    error('Failed to open TIF file.');
            end;
        end;
        
        jobMemory(1, 2) = ceil(1.2 * 14 * unitX);
    end;
    
    currentTime = clock;
    timeString = [...
        num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
        '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
    parameterDatabase = [pwd filesep 'jobParameters.calculateDelta.' timeString '.mat'];
    
    save(parameterDatabase,...
        'timepoints', 'timeClusters', 'dffSampling', 'subOffset', 'meanFraction', 'scaling', 'forceZero', ...
        'referenceOffset', 'medianFlag', 'kernelSize', ...
        'inputString', 'header', 'footer', 'dataType', 'specimen', 'cameras', 'channels', ...
        'inputType', 'outputType', 'jobMemory');
    
    if localRun(1) ~= 1
        disp(' ');
        disp(['Estimated memory consumption per workstation core is ' num2str(jobMemory(2)) ' GB.']);
        disp(' ');
        if jobMemory(2) <= coreMemory && nTicks > 1
            disp(['Submitting parametric cluster job for ' num2str(nTicks) ' time cluster(s).']);
        else
            disp(['Submitting individual cluster job(s) for ' num2str(nTicks) ' time cluster(s).']);
        end;
        disp(' ');
        
        if jobMemory(2) <= coreMemory && nTicks > 1
            cmdFunction = ['calculateDelta(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
            cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                '/parametric:1-' num2str(nTicks) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
            [status, systemOutput] = system(cmd);
            disp(['System response: ' systemOutput]);
        else
            for t = 1:nTicks
                cmdFunction = ['calculateDelta(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
                cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                    '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
                [status, systemOutput] = system(cmd);
                disp(['Submitting time cluster ' num2str(t, '%.4d') ': ' systemOutput]);
            end;
        end;
    else
        disp(' ');
        
        if localRun(2) > 1 && nTicks > 1
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            matlabpool(localRun(2));
            
            disp(' ');
            parfor t = 1:nTicks
                calculateDelta(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            disp(' ');
        else
            for t = 1:nTicks
                calculateDelta(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
        end;
    end;
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;