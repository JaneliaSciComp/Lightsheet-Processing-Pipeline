%% parameters

timepoints   = 0:9000;
references   = 4490:4510;
dffSampling  = 5;      % time interval in which new dF/F reference stacks are computed
kernelSize   = 3;      % kernel size of median filter used in dF/F reference stack computation

inputString  = 'X:\SiMView2\14-01-04\Dme_L1_57C10-GCaMP641_0_20140104_114246.corrected\Results\TimeFused';
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

medianFlags = mod(0:(numel(timepoints) - 1), dffSampling) == 0;
existingResultsDetected = 0;

for currentTP = length(timepoints):-1:1
    if dataType == 0
        outputPath = [inputString '.registered\SPM' num2str(specimen, '%.2d') '\TM' num2str(timepoints(currentTP), '%.6d')];
    else
        outputPath = [inputString '.registered\' header '.TM' num2str(timepoints(currentTP), '%.6d') footer];
    end;
    fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString '.shifts.mat'];
    fullFilePath = [outputPath '\' fileName];
    
    if exist(fullFilePath, 'file') == 2
        timepoints(currentTP) = [];
        medianFlags(currentTP) = [];
        existingResultsDetected = 1;
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    if ~existingResultsDetected
        nReferences = numel(references);
        
        disp(' ');
        disp(['Generating average reference stack from ' num2str(nReferences) ' time points.']);
        
        outputPath = [inputString '.registered'];
        if exist(outputPath, 'dir') == 7
            warning('Default output folder already exists. Abort processing if existing data should not be overwritten.');
        else
            mkdir(outputPath);
        end;
        
        for t = references
            disp(['Reading data for time point ' num2str(t)]);
            
            if dataType == 0
                inputPath = [inputString '\SPM' num2str(specimen, '%.2d') '\TM' num2str(t, '%.6d')];
                fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(t, '%.6d') configurationString inputExtension];
            else
                inputPath = [inputString '\' header '.TM' num2str(t, '%.6d') footer];
                fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(t, '%.6d') configurationString '.fusedStack' inputExtension];
            end;
            fullFilePath = [inputPath '\' fileName];
            
            currentStack = single(readImage(fullFilePath));
            
            if t == references(1)
                referenceStack = currentStack;
                referenceProjectionXY = max(currentStack, [], 3);
                referenceProjectionXZ = squeeze(max(currentStack, [], 2));
                referenceProjectionYZ = squeeze(max(currentStack, [], 1));
            else
                referenceStack = referenceStack + currentStack;
                referenceProjectionXY = referenceProjectionXY + max(currentStack, [], 3);
                referenceProjectionXZ = referenceProjectionXZ + squeeze(max(currentStack, [], 2));
                referenceProjectionYZ = referenceProjectionYZ + squeeze(max(currentStack, [], 1));
            end;
        end;
        
        disp('Saving average reference stack to disk.');
        
        referencePath = [outputPath '\Reference'];
        if exist(referencePath, 'dir') ~= 7
            mkdir(referencePath);
        end;
        
        stackName = ['SPM' num2str(specimen, '%.2d') configurationString '.referenceStack' outputExtension];
        xyProjectionName = ['SPM' num2str(specimen, '%.2d') configurationString '.reference_xyProjection' outputExtension];
        xzProjectionName = ['SPM' num2str(specimen, '%.2d') configurationString '.reference_xzProjection' outputExtension];
        yzProjectionName = ['SPM' num2str(specimen, '%.2d') configurationString '.reference_yzProjection' outputExtension];
        fullStackPath = [referencePath '\' stackName];
        fullXYProjectionPath = [referencePath '\' xyProjectionName];
        fullXZProjectionPath = [referencePath '\' xzProjectionName];
        fullYZProjectionPath = [referencePath '\' yzProjectionName];
        
        referenceStack = uint16(referenceStack / nReferences);
        referenceProjectionXY = uint16(referenceProjectionXY / nReferences);
        referenceProjectionXZ = uint16(referenceProjectionXZ / nReferences);
        referenceProjectionYZ = uint16(referenceProjectionYZ / nReferences);
        
        writeImage(referenceStack, fullStackPath);
        writeImage(referenceProjectionXY, fullXYProjectionPath);
        writeImage(referenceProjectionXZ, fullXZProjectionPath);
        writeImage(referenceProjectionYZ, fullYZProjectionPath);
    else
        disp(' ');
        disp('Existing processing results detected. Skipping generation of average reference stack.');
    end;
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        if dataType == 0
            inputPath = [inputString '\SPM' num2str(specimen, '%.2d') '\TM' num2str(timepoints(1), '%.6d')];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString inputExtension];
        else
            inputPath = [inputString '\' header '.TM' num2str(timepoints(1), '%.6d') footer];
            fileName = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString '.fusedStack' inputExtension];
        end;
        fullFilePath = [inputPath '\' fileName];
        
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
        
        jobMemory(1, 2) = ceil(1.2 * 2 * unitX);
    end;
    
    currentTime = clock;
    timeString = [...
        num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
        '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
    parameterDatabase = [pwd '\jobParameters.registerStacks.' timeString '.mat'];
    
    save(parameterDatabase,...
        'timepoints', 'references', 'medianFlags', 'dffSampling', 'kernelSize', ...
        'inputString', 'header', 'footer', 'dataType', 'specimen', 'cameras', 'channels', ...
        'inputType', 'outputType', 'jobMemory');
    
    if localRun(1) ~= 1
        disp(' ');
        disp(['Estimated memory consumption per workstation core is ' num2str(jobMemory(2)) ' GB.']);
        disp(' ');
        if jobMemory(2) <= coreMemory && nTimepoints > 1
            disp(['Submitting parametric cluster job for ' num2str(nTimepoints) ' time point(s).']);
        else
            disp(['Submitting individual cluster job(s) for ' num2str(nTimepoints) ' time point(s).']);
        end;
        disp(' ');
        
        if jobMemory(2) <= coreMemory && nTimepoints > 1
            cmdFunction = ['registerStacks(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
            cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
            [status, systemOutput] = system(cmd);
            disp(['System response: ' systemOutput]);
        else
            for t = 1:nTimepoints
                cmdFunction = ['registerStacks(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
                cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                    '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
                [status, systemOutput] = system(cmd);
                disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
            end;
        end;
    else
        disp(' ');
        
        if localRun(2) > 1 && nTimepoints > 1
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            matlabpool(localRun(2));
            
            disp(' ');
            parfor t = 1:nTimepoints
                registerStacks(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            disp(' ');
        else
            for t = 1:nTimepoints
                registerStacks(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
        end;
    end;
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;