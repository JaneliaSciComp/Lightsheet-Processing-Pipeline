function clusterCS_fn(filename, timepoints_per_node, job_number)

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
%% parameters
%% keller parameters

inputRoot        = input_parameters.inputRoot;%'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Results' filesep 'MultiFused';
outputRoot       = input_parameters.outputRoot;%'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Results' filesep 'MultiFused.Corrected';
headerPattern    = input_parameters.headerPattern;%'Mmu_E1_CAGTAG1.TM??????_multiFused_blending' filesep;
filePattern      = input_parameters.filePattern;%'SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';
configRoot       = input_parameters.configRoot;%'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Scripts' filesep 'SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';

timepoints       = input_parameters.timepoints;%0:570;
dataType         = input_parameters.dataType;%1;      % 0: process projections only, 1: process stacks and projections
percentile       = input_parameters.percentile;%[1 10]; % [percentileFlag, subsampling], percentileFlag - 0: calculate true minimum, otherwise: use as percentile for minimum calculation

inputType        = input_parameters.inputType;%0;      % 0: input data in KLB format
                                                       % 1: input data in JP2 format
                                                       % 2: input data in TIF format
outputType       = input_parameters.outputType;%0;     % 0: output data saved in KLB format
                                                       % 1: output data saved in JP2 format
                                                       % 2: output data saved in TIF format

% configuration of drift correction
correctDrift     = input_parameters.correctDrift;%1;
referenceTime    = input_parameters.referenceTime;%305;
referenceROI     = input_parameters.referenceROI;%[];     % [xStart xStop; yStart yStop; zStart zStop], provide empty vector to enforce use of maximum dimensions

% configuration of intensity normalization
correctIntensity = input_parameters.correctIntensity;%1;

maxStampDigits   = input_parameters.maxStampDigits;%6;

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
%% initialization

if localRun(1) == 1 && localRun(2) == 0
    localRun(2) = feature('numcores');
    disp(' ');
    disp([num2str(localRun(2)) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

nTimepoints = numel(timepoints);
inputDatabase = cell(nTimepoints, 1);
outputDatabase = cell(nTimepoints, 2);

for n = 1:nTimepoints
    inputDatabase{n} = [inputRoot filesep '' headerPattern filePattern];
    outputDatabase{n, 1} = [outputRoot filesep '' headerPattern];
    outputDatabase{n, 2} = [outputRoot filesep '' headerPattern filePattern];
    for i = maxStampDigits:-1:1
        precision = ['%.' num2str(i) 'd'];
        
        positions1 = strfind(inputDatabase{n}, repmat('?', [1 i]));
        for k = 1:length(positions1)
            inputDatabase{n}(positions1(k):(positions1(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
        
        positions2 = strfind(outputDatabase{n, 1}, repmat('?', [1 i]));
        for k = 1:length(positions2)
            outputDatabase{n, 1}(positions2(k):(positions2(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
        
        positions3 = strfind(outputDatabase{n, 2}, repmat('?', [1 i]));
        for k = 1:length(positions3)
            outputDatabase{n, 2}(positions3(k):(positions3(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
    end;
end;

load([configRoot filesep 'dimensions.mat']);
load([configRoot filesep 'dimensionsMax.mat']);
load([configRoot filesep 'dimensionsDeltas.mat']);

if correctDrift
    load([configRoot filesep 'driftTable.mat']);
    
    xOffsets = driftTable(:, 2);
    yOffsets = driftTable(:, 3);
    zOffsets = driftTable(:, 4);
    
    referenceIndex = find(dimensions(:, 1) == referenceTime, 1);
    
    xOffsets = xOffsets(referenceIndex) - xOffsets;
    yOffsets = yOffsets(referenceIndex) - yOffsets;
    zOffsets = zOffsets(referenceIndex) - zOffsets;
else
    xOffsets = [];
    yOffsets = [];
    zOffsets = [];
end;

if correctIntensity
    load([configRoot filesep 'intensityBackgrounds.mat']);
    load([configRoot filesep 'intensityFactors.mat']);
else
    intensityBackgrounds = [];
    intensityFactors = [];
end;

%% job submission

for t = numel(timepoints):-1:1
    if exist([outputDatabase{t, 2} '_yzProjection.corrected' outputExtension], 'file') == 2
        timepoints(t) = [];
        inputDatabase(t) = [];
        outputDatabase(t, :) = [];
    end;
end;

if ~isempty(timepoints)
    nTimepoints = numel(timepoints);
    
    if jobMemory(1) == 1 && localRun(1) ~= 1
        if dataType == 0
            unitX = 0;
        else
            unitX = 2 * dimensionsMax(1) * dimensionsMax(2) * dimensionsMax(3) / (1024 ^ 3);
        end;
        jobMemory(1, 2) = ceil(1.2 * 3 * unitX);
    end;
    
    try
        % % %     if localRun(1) ~= 1
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
        % % %         if jobMemory(2) <= coreMemory && nTimepoints > 1
        % % %             cmdFunction = ['correctStack(''' parameterDatabase ''', *, ' num2str(jobMemory(2)) ')'];
        % % %             cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
        % % %                 '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
        % % %             [status, systemOutput] = system(cmd);
        % % %             disp(['System response: ' systemOutput]);
        % % %         else
        % % %             for t = 1:nTimepoints
        % % %                 cmdFunction = ['correctStack(''' parameterDatabase ''', ' num2str(t) ', ' num2str(jobMemory(2)) ')'];
        % % %                 cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
        % % %                     '/progressmsg:"' num2str(jobMemory(2)) '" runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
        % % %                 [status, systemOutput] = system(cmd);
        % % %                 disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ': ' systemOutput]);
        % % %             end;
        % % %         end;
        % % %     else
        disp(' ');
        
        if localRun(2) > 1 && nTimepoints > 1
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            matlabpool(localRun(2));
            
            disp(' ');
            
            parfor t = 1:nTimepoints
                disp(['Submitting time point ' num2str(timepoints(t), '%.4d') ' to a local worker.']);
                correctStack(inputDatabase, outputDatabase, dataType, timepoints, dimensions, dimensionsMax, dimensionsDeltas, ...
                    correctIntensity, intensityBackgrounds, intensityFactors, correctDrift, inputType, outputType, ...
                    percentile, referenceROI, xOffsets, yOffsets, zOffsets, jobMemory, t, jobMemory(2));
            end;
            
            disp(' ');
            
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            
            disp(' ');
        else
            for t = 1:nTimepoints
                disp(['Processing time point ' num2str(timepoints(t), '%.4d')]);
                correctStack(inputDatabase, outputDatabase, dataType, timepoints, dimensions, dimensionsMax, dimensionsDeltas, ...
                    correctIntensity, intensityBackgrounds, intensityFactors, correctDrift, inputType, outputType, ...
                    percentile, referenceROI, xOffsets, yOffsets, zOffsets, jobMemory, t, jobMemory(2));
                disp(' ');
            end;
        end;
        % % %     end;
    catch ME
        rethrow(ME);
    end
else
    disp(' ');
    disp('Existing processing results detected for all selected time points. Submission aborted.');
    disp(' ');
end;