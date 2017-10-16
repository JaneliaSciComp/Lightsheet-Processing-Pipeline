%% parameters

timepoints   = 0:350;

inputDir     = 'R:\SV1\KM_15-12-05\Mmu_E1_mKate2\Results\TimeFused.Corrected\';
outputDir    = 'R:\SV1\KM_15-12-05\Mmu_E1_mKate2\Results\TimeFused.Corrected.Masks\';
header       = 'Mmu_E1_mKate2';
footer       = '_timeFused_blending';
stackLabel   = '.corrected.shifted';

specimen     = 0;
cameras      = 0:1;
channels     = 0;

threshold    = 25;     % threshold for binarizing volume
removeDirt   = 1;      % flag for removing all but the largest connected component determined by thresholding
connectivity = 6;      % three-dimensional connectivty for object removal (6, 18, 26)

localRun     = [1 8];  % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                       % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                       %         note: use "0" to enable automated detection of available CPU cores

inputType    = 0;      % 0: input data in KLB format
                       % 1: input data in JP2 format
                       % 2: input data in TIF format
outputType   = 0;      % 0: output data saved in KLB format
                       % 1: output data saved in JP2 format
                       % 2: output data saved in TIF format

jobMemory    = [2 0];  % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
                       % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                       %         note 1: slot 2 is only evaluated if automated memory management is disabled
                       %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)

%% job submission

if localRun(1) == 1 && localRun(2) == 0
    localRun(2) = min(12, feature('numcores'));
    disp(' ');
    disp([num2str(localRun(2)) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

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

for currentTP = length(timepoints):-1:1
    outputPath = [outputDir header '.TM' num2str(timepoints(currentTP), '%.6d') footer];
    fileNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') configurationString];
    outputName = [outputPath '\' fileNameHeader '.fusedStack_yzProjection.mask' outputExtension];
    
    if exist(outputName, 'file') == 2
        timepoints(currentTP) = [];
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        fileNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') configurationString];
        fileName = [inputDir header '.TM' num2str(timepoints(1), '%.6d') footer '\' fileNameHeader '.fusedStack' stackLabel inputExtension];
        
        try
            switch inputType
                case 0
                    headerInformation = readKLBheader(fileName);
                    stackDimensions = headerInformation.xyzct(1:3);
                case 1
                    [stackDimensions, bitDepth] = readJP2header(fileName);
                case 2
                    headerInformation = imfinfo(fileName);
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
        
        jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
    end;
    
    if jobMemory(1) ~= 2 && localRun(1) ~= 1
        disp(' ');
        disp(['Estimated memory consumption per workstation core is ' num2str(jobMemory(2)) ' GB.']);
        disp(' ');
        if jobMemory(2) <= coreMemory && nTimepoints > 1
            disp(['Submitting parametric cluster job for ' num2str(nTimepoints) ' time point(s).']);
        else
            disp(['Submitting individual cluster job(s) for ' num2str(nTimepoints) ' time point(s).']);
        end;
        disp(' ');
        
        currentTime = clock;
        timeString = [...
            num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
            '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
        parameterDatabase = [pwd '\jobParameters.createMask.' timeString '.mat'];
        
        save(parameterDatabase,...
            'timepoints', 'inputDir', 'outputDir', 'header', 'footer', 'stackLabel', 'specimen', 'cameras', 'channels', ...
            'threshold', 'removeDirt', 'connectivity', 'inputType', 'outputType', 'jobMemory');
        
        if jobMemory(2) <= coreMemory && nTimepoints > 1
            cmdFunction = ['createMask(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
            cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
            [status, systemOutput] = system(cmd);
            disp(['System response: ' systemOutput]);
        else
            for t = 1:nTimepoints
                cmdFunction = ['createMask(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
                cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                    '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
                [status, systemOutput] = system(cmd);
                disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
            end;
        end;
    else
        currentTime = clock;
        timeString = [...
            num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
            '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
        parameterDatabase = [pwd '\jobParameters.createMask.' timeString '.mat'];
        
        save(parameterDatabase,...
            'timepoints', 'inputDir', 'outputDir', 'header', 'footer', 'stackLabel', 'specimen', 'cameras', 'channels', ...
            'threshold', 'removeDirt', 'connectivity', 'inputType', 'outputType', 'jobMemory');
        
        disp(' ');
        
        if localRun(1) ~= 1
            for t = 1:nTimepoints
                fileNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') configurationString];
                fileName = [inputDir header '.TM' num2str(timepoints(t), '%.6d') footer '\' fileNameHeader '.fusedStack' stackLabel inputExtension];
                
                try
                    switch inputType
                        case 0
                            headerInformation = readKLBheader(fileName);
                            stackDimensions = headerInformation.xyzct(1:3);
                        case 1
                            [stackDimensions, bitDepth] = readJP2header(fileName);
                        case 2
                            headerInformation = imfinfo(fileName);
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
                
                jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
                
                disp(['Estimated memory consumption per workstation core at time point ' num2str(timepoints(t), '%.4d') ' is ' num2str(jobMemory(2)) ' GB.']);
                
                cmdFunction = ['createMask(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
                cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                    '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
                [status, systemOutput] = system(cmd);
                disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
            end;
        else
            if localRun(2) > 1 && nTimepoints > 1
                if matlabpool('size') > 0
                    matlabpool('close');
                end;
                matlabpool(localRun(2));
                
                disp(' ');
                
                parfor t = 1:nTimepoints
                    createMask(parameterDatabase, t, jobMemory(2));
                end;
                
                disp(' ');
                
                if matlabpool('size') > 0
                    matlabpool('close');
                end;
                
                disp(' ');
            else
                for t = 1:nTimepoints
                    createMask(parameterDatabase, t, jobMemory(2));
                end;
                
                disp(' ');
            end;
        end;
    end;
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;