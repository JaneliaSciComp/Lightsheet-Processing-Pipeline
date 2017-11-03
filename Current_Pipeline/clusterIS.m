%% parameters

timepoints   = 0:500;

rootFolder   = 'U:' filesep 'SiMView2.Processing' filesep '12-06-12' filesep 'Dme_E1_His2ARFP_01_20120612_175853.corrected' filesep 'Results' filesep 'TimeFused';
inputHeader1 = 'Dme_E1_His2ARFP.TM';
inputHeader2 = '_timeFused_blending';
inputHeader3 = 'SPC0_TM';
inputFooter  = '_CM0_CM1_CHN00_CHN01.fusedStack';
timeStamp    = '%.4d';

interMode    = [1 2];  % slot 1: 0 for 1-dimensional, 1 for 3-dimensional interpolation
                       % slot 2: 0 for linear, 1 for spline, 2 for cubic
scaling      = 2.031 / (6.5 / 16);
splitting    = 10;

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

for currentTP = length(timepoints):-1:1
    fileName = [...
        rootFolder '.Interpolated' filesep ...
        inputHeader1 num2str(timepoints(currentTP), timeStamp) inputHeader2 filesep '' ...
        inputHeader3 num2str(timepoints(currentTP), timeStamp) inputFooter inputExtension];
    
    if exist(fileName, 'file') == 2
        timepoints(currentTP) = [];
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    currentTime = clock;
    timeString = [...
        num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
        '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
    parameterDatabase = [pwd filesep 'jobParameters.interpolateStack.' timeString '.mat'];
    
    save(parameterDatabase,...
        'timepoints', 'rootFolder', 'inputHeader1', 'inputHeader2', 'inputHeader3', 'inputFooter', 'timeStamp', ...
        'interMode', 'scaling', 'splitting', 'inputType', 'outputType', 'jobMemory');
    
    if localRun(1) ~= 1
        if jobMemory(1) == 1
            fileName = [...
                rootFolder filesep '' ...
                inputHeader1 num2str(timepoints(1), timeStamp) inputHeader2 filesep '' ...
                inputHeader3 num2str(timepoints(1), timeStamp) inputFooter inputExtension];
            
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
                        warning('Failed to open KLB file header.');
                    case 1
                        warning('Failed to open JP2 file header.');
                    case 2
                        warning('Failed to open TIF file header.');
                end;
                try
                    stack = readImage(fileName);
                    unitX = 2 * size(stack, 1) * size(stack, 2) * size(stack, 3) / (1024 ^ 3);
                    clear stack;
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
            end;
            
            jobMemory(1, 2) = ceil(1.2 * 24 * unitX); % factor yet to be determined
        end;
        
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
            cmdFunction = ['interpolateStack(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
            cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
                '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
            [status, systemOutput] = system(cmd);
            disp(['System response: ' systemOutput]);
        else
            for t = 1:nTimepoints
                cmdFunction = ['interpolateStack(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
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
                interpolateStack(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
            
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            
            disp(' ');
        else
            for t = 1:nTimepoints
                interpolateStack(parameterDatabase, t, jobMemory(2));
            end;
            
            disp(' ');
        end;
    end;
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;