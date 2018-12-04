function clusterMF_fn(filename, timepoints_per_node, job_number)

input_parameters = loadjson(fileread(filename));
if ~isfield(input_parameters, 'verbose'), input_parameters.verbose = true; end
if input_parameters.verbose
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end
input_parameters = convert_limits_to_values(input_parameters, {'timepoints','xOffsets', 'yOffsets'});
if nargin == 3
    timepoints_per_node = str2double(timepoints_per_node);
    job_number = str2double(job_number);
    start_timepoint_index = timepoints_per_node*(job_number-1) + 1;
    end_timepoint_index = min(timepoints_per_node*job_number, numel(input_parameters.timepoints));
    input_parameters.timepoints = input_parameters.timepoints(start_timepoint_index:end_timepoint_index);
end
input_parameters.cropping = num2cell(input_parameters.cropping,size(input_parameters.cropping));
%% parameters

%% keller-provided parameters
% % % %% parameters

timepoints   = input_parameters.timepoints; 
inputString  = input_parameters.inputString;%'W:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected';
outputString = input_parameters.outputString;%'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'MultiFused' filesep 'Pha_E1_H2bRFP';
outputID     = input_parameters.outputID; %'_blending';
dataType     = input_parameters.dataType; %1;         % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks

specimen     = input_parameters.specimen;% 0;
angle        = input_parameters.angle;
cameras      = input_parameters.cameras;%0:1;
channels     = input_parameters.channels;%0:1;

reducedIO    = input_parameters.reducedIO; %1;         % 0 for full logging, 1 for minimal logging, 2 for slice logging
inputType    = input_parameters.inputType; %0;         % 0: input data in KLB format
                          % 1: input data in JP2 format
                          % 2: input data in TIF format
outputType   = input_parameters.outputType; %0;         % 0: output data saved in KLB format
                          % 1: output data saved in JP2 format
                          % 2: output data saved in TIF format
splitting    = input_parameters.splitting;
kernelSize   = input_parameters.kernelSize;
kernelSigma  = input_parameters.kernelSigma;
fraction     = input_parameters.fraction;  % minimal object size in binary mask, set to 0 to disable bwareaopen (default value 10 ^ -5)
maskMinimum  = input_parameters.maskMinimum;   % slot 1 provides percentile for minimum calculation (0 = true minimum), slot 2 provides subsampling for percentile computation
maskFactor   = input_parameters.maskFactor;       % fraction of (mean - minimum) intensity that is used to create the thresholded 3D mask
maskFusion   = input_parameters.maskFusion;         % 0 use only overlap regions for average mask, 1 combine full information from both contributing masks
padding      = input_parameters.padding;    % slot 1: set to 0 to disable padding of segmentation mask, 1 to replace zeros with center coordinate, 2 for smooth padding towards center coordinate
                          % slot 2: radius of disk element used to define size of transition zone in smooth padding mode (2)
slabSizes    = input_parameters.slabSizes;     % adaptive slab size for correlation, set to 0 to use slice mode (first entry: channels, second entry: cameras)
intSizes     = input_parameters.intSizes;    % adaptive slab size for intensity correction, set to 0 to use slice mode (first entry: channels, second entry: cameras)
percentile   = input_parameters.percentile;         % histogram position for determining background levels
subSampling  = input_parameters.subSampling;   % data set subsampling for percentile computation (first entry: slices, second entry: stacks)
medianFilter = input_parameters.medianFilter;       % median filter range, set to 0 to disable
preciseGauss = input_parameters.preciseGauss;         % 0 for uint16 precision (with rounding errors), 1 for double precision (for low-intensity data sets)
gaussFilter  = input_parameters.gaussFilter;     % gauss filter kernel size and sigma, set to 0 to disable
optimizer    = input_parameters.optimizer;         % 1 for Newton method, 2 for Nelder-Mead simplex method, 3 for customized gradient descent (by Fernando Amat)
correction   = input_parameters.correction; % slots 1 and 2 provide intensity correction flags for multi-channel and multi-camera: 0 to disable, 1 for calculated values, 2 for manual values (slots 3 and 4)
                          % slots 3 and 4 provide manual correction factors for the transformed channel/camera with respect to the reference channel/camera (will be inversely applied if < 1)
fusionType   = input_parameters.fusionType;         % 0 for adaptive blending, 1 for geometrical blending, 2 for wavelet fusion, 3 for averaging
transitions  = input_parameters.transitions;     % only used for geometrical blending, slot 1 provides y-center, slot 2 provides z-center:
                          % 0 for geometrical center of stack, or plane index to define custom transition plane
blending     = input_parameters.blending;    % blending ranges (first entry: channels, second entry: cameras)

enforceFlag  = input_parameters.enforceFlag;
verbose      = input_parameters.verbose;         % 0 for minimal text output, 1 for timing text output, 2 for full text ouput

cropping     = input_parameters.cropping; % ImageJ xy-ROI, start plane, stop plane (camera 0/1); set all to 0 to disable cropping
scaling      = input_parameters.scaling;             % axial step size <divided by> (pixel pitch <divided by> magnification)

leftFlags    = input_parameters.leftFlags;     % indicates in which channel the light sheet comes from the left-hand side (camera 0/1), 1 = reference channel, 2 = transformed channel

flipHFlag    = input_parameters.flipHFlag;         % indicates whether the stack recorded with the second camera should be flipped horizontally
flipVFlag    = input_parameters.flipVFlag;         % indicates whether the stack recorded with the second camera should be flipped vertically

xOffsets     = input_parameters.xOffsets; % x-shift of second data slice with respect to first data slice (in pixels)
yOffsets     = input_parameters.yOffsets; % y-shift of second data slice with respect to first data slice (in pixels)

frontFlag    = input_parameters.frontFlag;         % indicates for which camera the high-quality information is in the front, 1 = reference camera, 2 = transformed camera

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
%% job submission

if dataType == 1 && padding(1) > 0 
    warning('The use of padded clusterMF segmentation masks is inconsistent with the assignment of a segmented input data type.');
end;

if localRun(1) == 1 && localRun(2) == 0
    localRun(2) = feature('numcores');
    disp(' ');
    disp([num2str(localRun(2)) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

if length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
    processingMode = 0;
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
    processingMode = 1;
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
    processingMode = 2;
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion).');
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
    if exist([outputString '.TM' num2str(timepoints(currentTP), '%.6d') '_multiFused' outputID filesep...
            '/SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(currentTP), '%.6d') '_ANG' num2str(angle, '%.3d') configurationString '_jobCompleted.txt'], 'file') == 2
        timepoints(currentTP) = [];
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        inputFolder = [inputString filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(1), '%.6d')];
        header = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(1), '%.6d') '_ANG' num2str(angle, '%.3d')...
            '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d')];
        fileName = [inputFolder filesep '' header inputExtension];
        
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
        
        switch processingMode
            case 0
                if fusionType == 2
                    jobMemory(1, 2) = ceil(1.2 * 11 * unitX);
                else
                    jobMemory(1, 2) = ceil(1.2 * 8 * unitX);
                end;
            case 1
                if fusionType == 2
                    jobMemory(1, 2) = ceil(1.2 * 9 * unitX);
                else
                    jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
                end;
            case 2
                if fusionType == 2
                    switch dataType
                        case 0
                            jobMemory(1, 2) = ceil(1.2 * 9 * unitX);
                        case 1
                            jobMemory(1, 2) = ceil(1.2 * 7 * unitX);
                    end;
                else
                    switch dataType
                        case 0
                            jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
                        case 1
                            jobMemory(1, 2) = ceil(1.2 * 4 * unitX);
                    end;
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
% % %          uid = [char(java.net.InetAddress.getLocalHost.getHostName) '_' num2str(feature('getpid'))];
% % %         parameterDatabase = [pwd filesep 'jobParameters.multiFuse.' uid '_' num2str(input_parameters.timepoints(1))  '.mat'];
% % %         
% % %         save(parameterDatabase, ...
% % %             'timepoints', 'inputString', 'outputString', 'outputID', 'dataType', 'specimen', 'cameras', 'channels', 'reducedIO', 'inputType', 'outputType', ...
% % %             'splitting', 'kernelSize', 'kernelSigma', 'fraction', 'maskMinimum', 'maskFactor', 'maskFusion', 'padding', ...
% % %             'slabSizes', 'intSizes', 'percentile', 'subSampling', 'medianFilter', 'preciseGauss', 'gaussFilter', ...
% % %             'optimizer', 'correction', 'fusionType', 'transitions', 'blending', 'enforceFlag', 'verbose', ...
% % %             'cropping', 'scaling', 'leftFlags', 'flipHFlag', 'flipVFlag', 'xOffsets', 'yOffsets', 'frontFlag', 'jobMemory');
% % %         try
% % %             if jobMemory(2) <= coreMemory && nTimepoints > 1
% % %                 cmdFunction = ['multiFuse(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
% % %                 cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                     '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                 [status, systemOutput] = system(cmd);
% % %                 disp(['System response: ' systemOutput]);
% % %             else
% % %                 for t = 1:nTimepoints
% % %                     cmdFunction = ['multiFuse(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
% % %                     cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                         '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                     [status, systemOutput] = system(cmd);
% % %                     disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
% % %                 end;
% % %             end;
% % %         catch ME
% % %             rethrow(ME);
% % %         end
    else
        try
            if localRun(1) ~= 1
% % %                 uid = [char(java.net.InetAddress.getLocalHost.getHostName) '_' num2str(feature('getpid'))];
% % %                 parameterDatabase = [pwd filesep 'jobParameters.multiFuse.' uid '_' num2str(input_parameters.timepoints(1)) '.mat'];
% % %                 
% % %                 save(parameterDatabase, ...
% % %                     'timepoints', 'inputString', 'outputString', 'outputID', 'dataType', 'specimen', 'cameras', 'channels', 'reducedIO', 'inputType', 'outputType', ...
% % %                     'splitting', 'kernelSize', 'kernelSigma', 'fraction', 'maskMinimum', 'maskFactor', 'maskFusion', 'padding', ...
% % %                     'slabSizes', 'intSizes', 'percentile', 'subSampling', 'medianFilter', 'preciseGauss', 'gaussFilter', ...
% % %                     'optimizer', 'correction', 'fusionType', 'transitions', 'blending', 'enforceFlag', 'verbose', ...
% % %                     'cropping', 'scaling', 'leftFlags', 'flipHFlag', 'flipVFlag', 'xOffsets', 'yOffsets', 'frontFlag', 'jobMemory');
% % %                 
% % %                 disp(' ');
% % %                 for t = 1:nTimepoints
% % %                     inputFolder = [inputString filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d')];
% % %                     header = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') ...
% % %                         '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d')];
% % %                     fileName = [inputFolder filesep '' header inputExtension];
% % %                     
% % %                     try
% % %                         switch inputType
% % %                             case 0
% % %                                 headerInformation = readKLBheader(fileName);
% % %                                 stackDimensions = headerInformation.xyzct(1:3);
% % %                             case 1
% % %                                 [stackDimensions, bitDepth] = readJP2header(fileName);
% % %                             case 2
% % %                                 headerInformation = imfinfo(fileName);
% % %                                 stackDimensions = [headerInformation(1).Height headerInformation(1).Width numel(headerInformation)];
% % %                         end;
% % %                         unitX = 2 * stackDimensions(1) * stackDimensions(2) * stackDimensions(3) / (1024 ^ 3);
% % %                     catch errorMessage
% % %                         switch inputType
% % %                             case 0
% % %                                 error('Failed to open KLB file.');
% % %                             case 1
% % %                                 error('Failed to open JP2 file.');
% % %                             case 2
% % %                                 error('Failed to open TIF file.');
% % %                         end;
% % %                     end;
% % %                     
% % %                     switch processingMode
% % %                         case 0
% % %                             if fusionType == 2
% % %                                 jobMemory(1, 2) = ceil(1.2 * 11 * unitX);
% % %                             else
% % %                                 jobMemory(1, 2) = ceil(1.2 * 8 * unitX);
% % %                             end;
% % %                         case 1
% % %                             if fusionType == 2
% % %                                 jobMemory(1, 2) = ceil(1.2 * 9 * unitX);
% % %                             else
% % %                                 jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
% % %                             end;
% % %                         case 2
% % %                             if fusionType == 2
% % %                                 switch dataType
% % %                                     case 0
% % %                                         jobMemory(1, 2) = ceil(1.2 * 9 * unitX);
% % %                                     case 1
% % %                                         jobMemory(1, 2) = ceil(1.2 * 7 * unitX);
% % %                                 end;
% % %                             else
% % %                                 switch dataType
% % %                                     case 0
% % %                                         jobMemory(1, 2) = ceil(1.2 * 6 * unitX);
% % %                                     case 1
% % %                                         jobMemory(1, 2) = ceil(1.2 * 4 * unitX);
% % %                                 end;
% % %                             end;
% % %                     end;
% % %                     
% % %                     disp(['Estimated memory consumption per workstation core at time point ' num2str(timepoints(t), '%.4d') ' is ' num2str(jobMemory(2)) ' GB.']);
% % %                     
% % %                     cmdFunction = ['multiFuse(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
% % %                     cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
% % %                         '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
% % %                     [status, systemOutput] = system(cmd);
% % %                     disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
% % %                 end;
            else
                
                if localRun(2) > 1 && nTimepoints > 1
                    if matlabpool('size') > 0
                        matlabpool('close');
                    end;
                    matlabpool(localRun(2));
                    
                    disp(' ');
                    
                    parfor t = 1:nTimepoints
                        multiFuse(timepoints, inputString, outputString, outputID, dataType, specimen, angle, cameras, channels, reducedIO, inputType, outputType, ...
                            splitting, kernelSize, kernelSigma, fraction, maskMinimum, maskFactor, maskFusion, padding, ...
                            slabSizes, intSizes, percentile, subSampling, medianFilter, preciseGauss, gaussFilter, ...
                            optimizer, correction, fusionType, transitions, blending, enforceFlag, verbose, ...
                            cropping, scaling, leftFlags, flipHFlag, flipVFlag, xOffsets, yOffsets, frontFlag, jobMemory, t, jobMemory(2));
                    end;
                    
                    disp(' ');
                    
                    if matlabpool('size') > 0
                        matlabpool('close');
                    end;
                    
                    disp(' ');
                else
                    for t = 1:nTimepoints
                        multiFuse(timepoints, inputString, outputString, outputID, dataType, specimen, angle, cameras, channels, reducedIO, inputType, outputType, ...
                            splitting, kernelSize, kernelSigma, fraction, maskMinimum, maskFactor, maskFusion, padding, ...
                            slabSizes, intSizes, percentile, subSampling, medianFilter, preciseGauss, gaussFilter, ...
                            optimizer, correction, fusionType, transitions, blending, enforceFlag, verbose, ...
                            cropping, scaling, leftFlags, flipHFlag, flipVFlag, xOffsets, yOffsets, frontFlag, jobMemory, t, jobMemory(2));
                    end;
                    
                    disp(' ');
                end;
            end;
        catch ME
            rethrow(ME);
        end
    end;
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;
end