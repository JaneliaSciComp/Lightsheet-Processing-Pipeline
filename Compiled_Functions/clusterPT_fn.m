function clusterPT_fn(filename, timepoints_per_node, job_number)
%% parameters
% read json input
%distcomp.feature( 'LocalUseMpiexec', false );
input_parameters = loadjson(fileread(filename));
if ~isfield(input_parameters, 'verbose'), input_parameters.verbose = true; end
if input_parameters.verbose
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end
input_parameters = convert_limits_to_values(input_parameters, {'timepoints'});
if nargin == 3
    timepoints_per_node = str2double(timepoints_per_node);
    job_number = str2double(job_number);
    start_timepoint_index = timepoints_per_node*(job_number-1) + 1;
    end_timepoint_index = min(timepoints_per_node*job_number, numel(input_parameters.timepoints));
    input_parameters.timepoints = input_parameters.timepoints(start_timepoint_index:end_timepoint_index);
end
% root folder with raw data
inputFolder  = input_parameters.inputFolder;
userDefinedOutputFolder='';
if isfield(input_parameters,'outputFolder'), userDefinedOutputFolder = input_parameters.outputFolder; end%DGA added this so that the output folder can be different from the input
outputLabel  = input_parameters.outputLabel;
specimen     = input_parameters.specimen;    % specimen index to be processed
angle        = input_parameters.angle;       % angle to be processed
timepoints   = input_parameters.timepoints;  % time points to be processed
cameras      = input_parameters.cameras;     % camera indices to be processed
channels     = input_parameters.channels;    % channel indices to be processed

dimensions   = input_parameters.dimensions;  % override parameter for user-provided raw stack dimensions (ImageJ-x, ImageJ-y, z/t)
                                             % note: assign empty vector [] to indicate that override is not needed (requires XML meta data)

% use highly saturated global XY projection to determine the following four parameters (can be left empty to enforce maximum ROI)
startsLeft   = input_parameters.startsLeft;  % cropping left start coordinates for each camera (ImageJ-x convention)
startsTop    = input_parameters.startsTop;   % cropping top start coordinates for each camera (ImageJ-y convention)
widths       = input_parameters.widths;      % cropping width for each camera (ImageJ-w convention)
heights      = input_parameters.heights;     % cropping height for each camera (ImageJ-h convention)

% use highly saturated global XZ projection to determine the following two parameters (can be left empty to enforce maximum ROI)
startsFront  = input_parameters.startsFront; % cropping front start coordinates for each camera (ImageJ-y convention)
depths       = input_parameters.depths;      % cropping depth for each camera (ImageJ-h convention)

inputType    = input_parameters.inputType;   % 0: input data in TIF format
                                             % 1: input data in JP2 format
                                             % 2: input data in binary stack format (normal experiment mode, structured)
                                             % 3: input data in binary stack format (normal experiment mode, unstructured)
                                             % 4: input data in binary stack format (HS single-plane experiment mode)
                                             %    note: pixel correction and segmentFlag are inactive in this mode
outputType   = input_parameters.outputType;  % 0: output data saved in KLB format
                                             % 1: output data saved in JP2 format
                                             % 2: output data saved in TIF format
correctTIFF  = input_parameters.correctTIFF; % 0: never transpose output data, 1: transpose output data if TIFF parity is odd

rotationFlag = input_parameters.rotationFlag;% 0: do not rotate image stacks, 1: rotate image stacks by 90 degrees clockwise, -1: rotate image stacks by 90 degrees counter-clockwise
medianRange  = input_parameters.medianRange; % kernel x/y-size for median filter for dead pixel detection and removal
percentile   = input_parameters.percentile;  % slot 1: background percentile for mask calculation
                                             % slot 2: background percentile for intensity correction
                                             % slot 3: volume sub-sampling for background estimation

segmentFlag  = input_parameters.segmentFlag; % 0: do not remove background in output stacks, 1: remove background in output stacks (default)
flipHFlag    = input_parameters.flipHFlag;   % indicates whether the stack recorded with the second camera should be flipped horizontally
flipVFlag    = input_parameters.flipVFlag;   % indicates whether the stack recorded with the second camera should be flipped vertically
splitting    = input_parameters.splitting;   % level of stack splitting when performing Gauss convolution
kernelSize   = input_parameters.kernelSize;  % Gauss kernel size
kernelSigma  = input_parameters.kernelSigma; % Gauss kernel sigma
scaling      = input_parameters.scaling;     % axial step size <divided by> (pixel pitch <divided by> magnification)
references   = input_parameters.references;  % list of reference channel groups that define segmentation masks for (content-)dependent channel groups
                                             % note 1: leave empty if all channels are to be treated as independent channels
                                             % note 2: provide multiple elements per row (i.e. a channel group) to fuse segmentation masks
dependents   = input_parameters.dependents;  % list of (content-)dependent channel groups that relate to each reference channel group
thresholds   = input_parameters.thresholds;  % adaptive thresholds for each reference channel group
                                             % note: provide single threshold if "references" is left empty

loggingFlag  = input_parameters.loggingFlag; % 0: do not save processing logs (default), 1: save processing logs
verbose      = input_parameters.verbose;     % 0: do not display processing information, 1: display processing information

% Deprecated, unused values that are just set here to defaults as
% placeholders until code is changed
localRun     = [1 0];                        % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                                             %         note: cluster submission requires a Windows cluster with HPC 2008 support (see "job submit" commands below)
                                             % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                                             %         note: use "0" to enable automated detection of available CPU cores

jobMemory    = [1 0];                        % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
                                             %         note: setting jobMemory(1) to "2" is incompatible with inputType == 4 or numel(dimensions) ~= 0
                                             % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                                             %         note 1: slot 2 is only evaluated if automated memory management is disabled
                                             %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)
                                                       % note: parameter is only required for cluster submission

%% keller-provided inputs
% % %inputFolder  = ['W:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802'];
% % %outputLabel  = '';
% % % specimen     = 0;                   % specimen index to be processed
% % % timepoints   = 0:733;               % time points to be processed
% % % cameras      = [0 1];               % camera indices to be processed
% % % channels     = [0 1];               % channel indices to be processed
% % % 
% % % dimensions   = [603 1272 125]; %[];                  % override parameter for user-provided raw stack dimensions (ImageJ-x, ImageJ-y, z/t)
% % %                                     % note: assign empty vector [] to indicate that override is not needed
% % % 
% % % % use highly saturated global XY projection to determine the following four parameters (can be left empty to enforce maximum ROI)
% % % startsLeft   = [ 372  216];         % cropping left start coordinates for each camera (ImageJ-x convention)
% % % startsTop    = [  68   55];         % cropping top start coordinates for each camera (ImageJ-y convention)
% % % widths       = [1148 1148];         % cropping width for each camera (ImageJ-w convention)
% % % heights      = [1468 1468];         % cropping height for each camera (ImageJ-h convention)
% % % 
% % % % use highly saturated global XZ projection to determine the following two parameters (can be left empty to enforce maximum ROI)
% % % startsFront  = [   0    0];         % cropping front start coordinates for each camera (ImageJ-y convention)
% % % depths       = [ 234  234];         % cropping depth for each camera (ImageJ-h convention)
% % % 
% % % inputType    = 0;                   % 0: input data in TIFF format
% % %                                     % 1: input data in JP2 format
% % %                                     % 2: input data in binary stack format (normal experiment mode, structured)
% % %                                     % 3: input data in binary stack format (normal experiment mode, unstructured)
% % %                                     % 4: input data in binary stack format (HS single-plane experiment mode)
% % %                                     %    note: pixel correction and segmentFlag are inactive in this mode
% % % outputType   = 0;                   % 0: output data saved in KLB format
% % %                                     % 1: output data saved in JP2 format
% % %                                     % 2: output data saved in TIF format
% % % correctTIFF  = 0;                   % 0: never transpose output data, 1: transpose output data if TIFF parity is odd
% % % 
% % % medianRange  = [3 3];               % kernel x/y-size for median filter for dead pixel detection and removal
% % %                                     % note: set to [] to disable dead pixel detection and removal
% % % percentile   = [1 5 100];           % slot 1: background percentile for mask calculation
% % %                                     % slot 2: background percentile for intensity correction
% % %                                     % slot 3: volume sub-sampling for background estimation
% % % 
% % % segmentFlag  = 1;                   % 0: do not remove background in output stacks
% % %                                     % 1: remove background in output stacks (default)
% % %                                     % 2: compute and save only 3D segmentation masks (first step of global masking workflow);
% % %                                     %    once all masks have been generated, rerunning clusterPT will generate a global mask (second step of global masking workflow)
% % %                                     % 3: apply global 3D segmentation mask to all image stacks (third step of global masking workflow)
% % % rotationFlag = 0;                   % 0: do not rotate image stacks, 1: rotate image stacks by 90 degrees clockwise, -1: rotate image stacks by 90 degrees counter-clockwise
% % % flipHFlag    = 0;                   % indicates whether the stack recorded with the second camera should be flipped horizontally
% % % flipVFlag    = 0;                   % indicates whether the stack recorded with the second camera should be flipped vertically
% % % splitting    = 10;                  % level of stack splitting when performing Gauss convolution
% % % kernelSize   = 5;                   % Gauss kernel size
% % % kernelSigma  = 2;                   % Gauss kernel sigma
% % % scaling      = 2.031 / (6.5 / 16);  % axial step size <divided by> (pixel pitch <divided by> magnification)
% % % references   = [0; 1];              % list of reference channel groups that define segmentation masks for (content-)dependent channel groups
% % %                                     % note 1: leave empty if all channels are to be treated as independent channels
% % %                                     % note 2: provide multiple elements per row (i.e. a channel group) to fuse segmentation masks
% % % dependents   = [];                  % list of (content-)dependent channel groups that relate to each reference channel group
% % % thresholds   = [0.4 0.4];           % adaptive thresholds for each reference channel group
% % %                                     % note: provide single threshold if "references" is left empty
% % % 
% % % loggingFlag  = 0;                   % 0: do not save processing logs (default), 1: save processing logs
% % % verbose      = 1;                   % 0: do not display processing information, 1: display processing information
% % % 
% % % localRun     = [1 0];               % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
% % %                                     % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
% % %                                     %         note: use "0" to enable automated detection of available CPU cores
% % % 
% % % jobMemory    = [1 0];               % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
% % %                                     %         note: setting jobMemory(1) to "2" is incompatible with inputType == 4 or numel(dimensions) ~= 0
% % %                                     % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
% % %                                     %         note 1: slot 2 is only evaluated if automated memory management is disabled
% % %                                     %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode
% % % 
% % % coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)

%% initialization

if inputType == 4 && segmentFlag > 0
    error('Error: clusterPT segmentation features cannot be applied to image data recorded in HS single-plane experiment mode.');
end;

if segmentFlag == 2
    timepointsBackup = timepoints;
end;

if localRun(1) == 1 && localRun(2) == 0
    localRun(2) = feature('numcores');
    disp(' ');
    disp([num2str(localRun(2)) ' CPU cores were detected and will be allocated for parallel processing.']);
    disp(' ');
end;

% validate cropping region
if ~(isempty(startsTop) || isempty(heights) || isempty(startsLeft) || isempty(widths) || isempty(startsFront) || isempty(depths))
    nCameras = numel(cameras);
    if nCameras > 1
        for i = 2:nCameras
            if widths(1) ~= widths(i) || heights(1) ~= heights(i) || depths(1) ~= depths(i)
                warning(['Cropping regions for camera ' num2str(cameras(1)) ' and camera ' num2str(cameras(i)) ' are not the same size.']);
            end;
        end;
    end;
else
    warning('Incomplete cropping information provided by user. Meta data from XML files will be used to determine maximum ROI size.');
end;

%DGA: Added the following so the user can define the output folder
if ~isempty(userDefinedOutputFolder)
    templateOutputFolder = userDefinedOutputFolder;
else
    templateOutputFolder = inputFolder;
end

% create task-specific input and output folders
if ~isempty(outputLabel)
    projectionFolder = [templateOutputFolder '.corrected.' outputLabel '.projections'];
else
    projectionFolder = [templateOutputFolder '.corrected.projections'];
end;

if inputType == 4
    if ~isempty(outputLabel)
        outputFolder = [templateOutputFolder '.corrected.' outputLabel];
    else
        outputFolder = [templateOutputFolder '.corrected'];
    end;
else
    if ~isempty(outputLabel)
        outputFolder = [templateOutputFolder '.corrected.' outputLabel filesep 'SPM' num2str(specimen, '%.2d')];
    else
        outputFolder = [templateOutputFolder '.corrected' filesep 'SPM' num2str(specimen, '%.2d')];
    end;
end;
if ~isempty(outputLabel)
    globalMaskFolder = [templateOutputFolder '.globalMask.' outputLabel];
else
    globalMaskFolder = [templateOutputFolder '.globalMask'];
end;
if inputType ~= 3 && inputType ~= 4
    inputFolder = [inputFolder filesep 'SPM' num2str(specimen, '%.2d')];
end;

if exist(outputFolder, 'dir') == 7 && exist(projectionFolder, 'dir') == 7
    warning('Default output folders already exist. Abort processing if existing data should not be overwritten.');
else
    mkdir(outputFolder);
    if inputType ~= 4 && exist(projectionFolder, 'dir') ~= 7
        mkdir(projectionFolder);
    end;
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

% evalute background files and copy to output folder
if inputType == 0 || inputType == 2 || inputType == 3 || inputType == 4
    backgroundFiles = dir([inputFolder filesep 'Background_*.tif']);
else % inputType == 1
    backgroundFiles = dir([inputFolder filesep 'Background_*.jp2']);
end;
if isempty(backgroundFiles)
    error('Error: Background files are missing.');
end;
backgroundValues = zeros(length(backgroundFiles), 1);
for i = 1:length(backgroundFiles)
    source = [inputFolder filesep '' backgroundFiles(i).name];
    
    backgroundImage = readImage(source);
    backgroundValues(i) = single(prctile(backgroundImage(:), 3));
    
    if rotationFlag == 1
        backgroundImage = rot90(backgroundImage, 3);
    elseif rotationFlag == -1
        backgroundImage = rot90(backgroundImage, 1);
    end;
    
    target = [outputFolder filesep '' backgroundFiles(i).name(1:(end - 4)) outputExtension];
    writeImage(backgroundImage, target);
end;

%% job submission

for currentTP = length(timepoints):-1:1
    if inputType == 4
        currentMissingFlag = 0;
        for c = cameras
            for h = channels
                % check for presence of final image file written by each processTimepoint.m process
                outputFile = [outputFolder filesep 'Sequence' num2str(timepoints(currentTP), '%.6d') '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') outputExtension];
                if exist(outputFile, 'file') ~= 2
                    currentMissingFlag = 1;
                    break;
                end;
            end;
        end;
    else
        currentHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') '_ANG' num2str(angle, '%.3d')];
        currentMissingFlag = 0;
        for c = cameras
            if segmentFlag == 2
                if ~isempty(references)
                    for h = references(1, :)
                        % check for presence of final image file written by each processTimepoint.m process
                        outputFile = [outputFolder filesep 'TM' num2str(timepoints(currentTP), '%.6d') filesep currentHeader ...
                            '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.segmentationMask' outputExtension];
                        if exist(outputFile, 'file') ~= 2
                            currentMissingFlag = 1;
                            break;
                        end;
                    end;
                else
                    for h = channels
                        % check for presence of final image file written by each processTimepoint.m process
                        outputFile = [outputFolder filesep 'TM' num2str(timepoints(currentTP), '%.6d') filesep currentHeader ...
                            '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.segmentationMask' outputExtension];
                        if exist(outputFile, 'file') ~= 2
                            currentMissingFlag = 1;
                            break;
                        end;
                    end;
                end;
            else
                for h = channels
                    % check for presence of final image file written by each processTimepoint.m process
                    outputFile = [projectionFolder filesep '' currentHeader '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.yzProjection' outputExtension];
                    if exist(outputFile, 'file') ~= 2
                        currentMissingFlag = 1;
                        break;
                    end;
                end;
            end;
        end;
    end;
    if currentMissingFlag == 0
        timepoints(currentTP) = [];
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        if numel(dimensions) == 3 && (dimensions(1) * dimensions(2) * dimensions(3) ~= 0)
            unitX = 2 * dimensions(1) * dimensions(2) * dimensions(3) / (1024 ^ 3);
        else
            if inputType == 3 || inputType == 4
                xmlName = [inputFolder filesep 'ch' num2str(channels(1)) '.xml'];
            else
                xmlName = [inputFolder filesep 'TM' num2str(timepoints(1), '%.5d') filesep 'ch' num2str(channels(1)) '.xml'];
            end;
            
            try
                stackDimensions = [];
                
                fid = fopen(xmlName, 'r');
                if fid < 1
                    error('Could not open file %s for reading.', xmlName);
                end;
                while true
                    s = fgetl(fid);
                    if ~ischar(s)
                        break;
                    end;
                    m = regexp(s, '<info dimensions="([^#]*)"', 'tokens', 'once');
                    if ~isempty(m)
                        stackDimensions = m{1};
                        break;
                    end;
                end;
                fclose(fid);
                
                if isempty(stackDimensions)
                    error('Unable to retrieve stack dimensions from file %s.', xmlName);
                end
                stackDimensions = str2double(regexp(stackDimensions, '[^\d]+', 'split'));
                stackDimensions = reshape(stackDimensions, [3, numel(stackDimensions) / 3])';
                if any(isnan(stackDimensions))
                    error('Unable to correctly parse stack dimensions retrieved from file %s.', xmlName);
                end;
                
                unitX = 2 * stackDimensions(1, 1) * stackDimensions(1, 2) * stackDimensions(1, 3) / (1024 ^ 3);
            catch errorMessage
                error('Failed to open XML file.');
            end;
        end;
        
        if inputType == 4 % for HS single-plane experiment mode, the estimated memory consumption is 1 * unitX (+10% safety)
            if rotationFlag == 0
                jobMemory(1, 2) = ceil(1.2 * unitX);
            elseif rotationFlag == -1
                jobMemory(1, 2) = ceil(1.2 * 2 * unitX);
            elseif rotationFlag == 1
                jobMemory(1, 2) = ceil(1.2 * 3 * unitX);
            end;
        else              % otherwise, the estimated memory consumption is (2*n + 2) * unitX (n = max. number of channels used in parallel, +10% safety)
            if isempty(references)
                nReferences = 1;
            else
                nReferences = 0;
                for r = 1:size(references, 1)
                    nReferences = max(nReferences, numel(references(r, :)));
                end;
            end;
            if rotationFlag == 0
                jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 2) * unitX);
            elseif rotationFlag == -1
                jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 3) * unitX);
            elseif rotationFlag == 1
                jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 4) * unitX);
            end;
        end;
    end;
    
    if jobMemory(1) ~= 2 && localRun(1) ~= 1
% % %         disp(' ');
% % %         disp(['Estimated memory consumption per workstation core is ' num2str(jobMemory(2)) ' GB.']);
% % %         disp(' ');
% % %         if jobMemory(2) <= coreMemory && nTimepoints > 1
% % %             disp(['Submitting parametric cluster job for ' num2str(nTimepoints) ' time point(s).']);
% % %         else
% % %             disp(['Submitting individual cluster job(s) for ' num2str(nTimepoints) ' time point(s).']);
% % %         end;
% % %         disp(' ');
% % %         
% % %         uid = [char(java.net.InetAddress.getLocalHost.getHostName) '_' num2str(feature('getpid'))];
% % %         parameterDatabase = [pwd filesep 'jobParameters.processTimepoint.' uid '_' num2str(input_parameters.timepoints(1))  '.mat'];
% % %         
% % %         save(parameterDatabase, ...
% % %             'inputFolder', 'outputFolder', 'projectionFolder', 'globalMaskFolder', 'specimen', 'timepoints', 'cameras', 'channels', 'dimensions', ...
% % %             'startsLeft', 'startsTop', 'widths', 'heights', 'startsFront', 'depths', 'inputType', 'outputType', 'correctTIFF', 'rotationFlag', ...
% % %             'medianRange', 'percentile', 'segmentFlag', 'flipHFlag', 'flipVFlag', 'splitting', 'kernelSize', 'kernelSigma', 'scaling', ...
% % %             'references', 'dependents', 'thresholds', 'loggingFlag', 'verbose', 'backgroundValues', 'jobMemory');
% % %         try
% % %             if jobMemory(2) <= coreMemory && nTimepoints > 1
% % %                 cmdFunction = ['processTimepoint(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
% % %                 cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                     '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                 [status, systemOutput] = system(cmd);
% % %                 disp(['System response: ' systemOutput]);
% % %             else
% % %                 for t = 1:nTimepoints
% % %                     cmdFunction = ['processTimepoint(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
% % %                     cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                         '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                     [status, systemOutput] = system(cmd);
% % %                     disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
% % %                 end;
% % %             end;
% % %         catch ME
% % %             error(ME.message);
% % %         end
    else
        if localRun(1) ~= 1
            if inputType == 4
                error('Inconsistent parameters: "jobMemory(1) = 2" and "inputType = 4".');
            end;
            if numel(dimensions) == 3 && (dimensions(1) * dimensions(2) * dimensions(3) ~= 0)
                error('Inconsistent parameters: "jobMemory(1) = 2" is in conflict with assignment of non-empty "dimensions" vector.');
            end;
        end;
        
        if isempty(references)
            nReferences = 1;
        else
            nReferences = 0;
            for r = 1:size(references, 1)
                nReferences = max(nReferences, numel(references(r, :)));
            end;
        end;
        
        
        disp(' ');
        try
            if localRun(1) ~= 1
% % %                 uid = [char(java.net.InetAddress.getLocalHost.getHostName) '_' num2str(feature('getpid'))];
% % %                 parameterDatabase = [pwd filesep 'jobParameters.processTimepoint.' uid '_' num2str(input_parameters.timepoints(1)) '.mat'];
% % %                 
% % %                 save(parameterDatabase, ...
% % %                     'inputFolder', 'outputFolder', 'projectionFolder', 'globalMaskFolder', 'specimen', 'timepoints', 'cameras', 'channels', 'dimensions', ...
% % %                     'startsLeft', 'startsTop', 'widths', 'heights', 'startsFront', 'depths', 'inputType', 'outputType', 'correctTIFF', 'rotationFlag', ...
% % %                     'medianRange', 'percentile', 'segmentFlag', 'flipHFlag', 'flipVFlag', 'splitting', 'kernelSize', 'kernelSigma', 'scaling', ...
% % %                     'references', 'dependents', 'thresholds', 'loggingFlag', 'verbose', 'backgroundValues', 'jobMemory');
% % %                 for t = 1:nTimepoints
% % %                     xmlName = [inputFolder filesep 'TM' num2str(timepoints(t), '%.5d') filesep 'ch' num2str(channels(1)) '.xml'];
% % %                     
% % %                     try
% % %                         stackDimensions = [];
% % %                         
% % %                         fid = fopen(xmlName, 'r');
% % %                         if fid < 1
% % %                             error('Could not open file %s for reading.', xmlName);
% % %                         end;
% % %                         while true
% % %                             s = fgetl(fid);
% % %                             if ~ischar(s)
% % %                                 break;
% % %                             end;
% % %                             m = regexp(s, '<info dimensions="([^#]*)"', 'tokens', 'once');
% % %                             if ~isempty(m)
% % %                                 stackDimensions = m{1};
% % %                                 break;
% % %                             end;
% % %                         end;
% % %                         fclose(fid);
% % %                         
% % %                         if isempty(stackDimensions)
% % %                             error('Unable to retrieve stack dimensions from file %s.', xmlName);
% % %                         end
% % %                         stackDimensions = str2double(regexp(stackDimensions, '[^\d]+', 'split'));
% % %                         stackDimensions = reshape(stackDimensions, [3, numel(stackDimensions) / 3])';
% % %                         if any(isnan(stackDimensions))
% % %                             error('Unable to correctly parse stack dimensions retrieved from file %s.', xmlName);
% % %                         end;
% % %                         
% % %                         unitX = 2 * stackDimensions(1, 1) * stackDimensions(1, 2) * stackDimensions(1, 3) / (1024 ^ 3);
% % %                     catch errorMessage
% % %                         error('Failed to open XML file.');
% % %                     end;
% % %                     
% % %                     % inputMode ~= 3, the estimated memory consumption is (2*n + 2) * unitX (n = max. number of channels used in parallel, +10% safety)
% % %                     if rotationFlag == 0
% % %                         jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 2) * unitX);
% % %                     elseif rotationFlag == -1
% % %                         jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 3) * unitX);
% % %                     elseif rotationFlag == 1
% % %                         jobMemory(1, 2) = ceil(1.2 * (2 * nReferences + 4) * unitX);
% % %                     end;
% % %                     
% % %                     disp(['Estimated memory consumption per workstation core at time point ' num2str(timepoints(t), '%.4d') ' is ' num2str(jobMemory(2)) ' GB.']);
% % %                     
% % %                     cmdFunction = ['processTimepoint(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
% % %                     cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                         '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                     [status, systemOutput] = system(cmd);
% % %                     disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
% % %                end;
            else
                if localRun(2) > 1 && nTimepoints > 1
                    if matlabpool('size') > 0
                        matlabpool('close');
                    end;
                    matlabpool(localRun(2));
                    
                    disp(' ');
                    
                    parfor t = 1:nTimepoints
                        disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ' to a local worker.']);
                        processTimepoint(inputFolder, outputFolder, projectionFolder, globalMaskFolder, specimen, angle, timepoints, cameras, channels, dimensions, ...
                            startsLeft, startsTop, widths, heights, startsFront, depths, inputType, outputType, correctTIFF, rotationFlag, ...
                            medianRange, percentile, segmentFlag, flipHFlag, flipVFlag, splitting, kernelSize, kernelSigma, scaling, ...
                            references, dependents, thresholds, loggingFlag, verbose, backgroundValues, jobMemory, t, jobMemory(2));
                    end;
                    
                    disp(' ');
                    
                    if matlabpool('size') > 0
                        matlabpool('close');
                    end;
                    
                    disp(' ');
                else
                    for t = 1:nTimepoints
                        disp(['Processing time point ' num2str(timepoints(t), '%.4d')]);
                        processTimepoint(inputFolder, outputFolder, projectionFolder, globalMaskFolder, specimen, angle, timepoints, cameras, channels, dimensions, ...
                            startsLeft, startsTop, widths, heights, startsFront, depths, inputType, outputType, correctTIFF, rotationFlag, ...
                            medianRange, percentile, segmentFlag, flipHFlag, flipVFlag, splitting, kernelSize, kernelSigma, scaling, ...
                            references, dependents, thresholds, loggingFlag, verbose, backgroundValues, jobMemory, t, jobMemory(2));
                        disp(' ');
                    end;
                end;
            end;
        catch ME
            rethrow(ME);
        end
    end;
elseif segmentFlag == 2
    disp(' ');
    disp('Existing processing results detected for all selected time points. Starting data crawler for construction of global segmentation mask.');
    
    if exist(globalMaskFolder, 'dir') ~= 7
        mkdir(globalMaskFolder);
    end;
    
    nTimepoints = numel(timepointsBackup);
    
    disp(' ');
    
    maskCount = 0;
    
    for t = 1:nTimepoints
        disp(['Reading masks for time point ' num2str(timepointsBackup(t), '%.4d')]);
        
        maskFolder = [outputFolder filesep 'TM' num2str(timepointsBackup(t), '%.6d') filesep ];
        maskHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepointsBackup(t), '%.6d') '_ANG' num2str(angle, '%.3d') '_CM*_CHN*.segmentationMask' outputExtension];
        maskFiles = dir([maskFolder maskHeader]);
        
        for i = 1:numel(maskFiles)
            disp(['* ' maskFiles(i).name]);
            if t == 1 && i == 1
                globalMask = readImage([maskFolder maskFiles(i).name]) == 1;
            else
                globalMask = globalMask | (readImage([maskFolder maskFiles(i).name]) == 1);
            end;
            maskCount = maskCount + 1;
        end;
    end;
    
    disp(' ');
    disp(['Total number of masks read: ' num2str(maskCount)]);
    disp('Writing global mask to disk');
    
    writeImage(uint16(globalMask), [globalMaskFolder filesep 'Global.segmentationMask' outputExtension]);
    
    xzSliceMask = zeros(size(globalMask, 1), size(globalMask, 3), 'uint16');
    for i = 1:splitting
        slabStart = round((i - 1) * size(globalMask, 1) / splitting + 1);
        slabStop  = round(i * size(globalMask, 1) / splitting);
        coordinateMask = repmat(1:size(globalMask, 2), [(slabStop - slabStart + 1) 1 size(globalMask, 3)]);
        coordinateMask(globalMask(slabStart:slabStop, :, :) == 0) = NaN;
        xzSliceMask(slabStart:slabStop, :) = uint16(round(squeeze(nanmean(coordinateMask, 2))));
    end;
    clear coordinateMask;
    
    xySliceMask = zeros(size(globalMask, 1), size(globalMask, 2), 'uint16');
    for i = 1:splitting
        slabStart = round((i - 1) * size(globalMask, 1) / splitting + 1);
        slabStop  = round(i * size(globalMask, 1) / splitting);
        coordinateMask = repmat(reshape(1:size(globalMask, 3), [1 1 size(globalMask, 3)]), [(slabStop - slabStart + 1) size(globalMask, 2) 1]);
        coordinateMask(globalMask(slabStart:slabStop, :, :) == 0) = NaN;
        xySliceMask(slabStart:slabStop, :) = uint16(round(squeeze(nanmean(coordinateMask, 3))));
    end;
    clear coordinateMask;
    
    writeImage(xzSliceMask, [globalMaskFolder filesep 'Global.xzMask' outputExtension]);
    writeImage(xySliceMask, [globalMaskFolder filesep 'Global.xyMask' outputExtension]);
    
    disp(' ');
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;
end