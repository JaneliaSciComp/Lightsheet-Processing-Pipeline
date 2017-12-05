% % % %% parameters

input_parameters.inputRoot        = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example\ Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.TimeFused/';
input_parameters.inputPattern     = ['Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.TM??????_timeFused_blending' filesep 'SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack'];
input_parameters.configRoot       = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Configurations/SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';

input_parameters.timepoints       = 0:50;
input_parameters.gamma            = 1;
input_parameters.percentile       = 1;

input_parameters.inputType        = 0;       % 0: input data in KLB format
                                             % 1: input data in JP2 format
                                             % 2: input data in TIF format
input_parameters.outputType       = 0;       % 0: output data saved in KLB format
                                             % 1: output data saved in JP2 format
                                             % 2: output data saved in TIF format

% configuration of intensity normalization
input_parameters.intensityFlag    = 1;       % flag for enabling/disabling intensity normalization
input_parameters.useStacks        = [1 10];  % 0: use projections for intensity estimate, 1: use stacks for intensity estimate (second parameter provides sub-sampling rate
input_parameters.histogramBins    = 0:(2^16 - 1);
input_parameters.threshold        = 10;      % threshold for histogram computation
input_parameters.backgroundSlot   = 0;       % 0: no background correction, 1: top right, 2: bottom right, 3: bottom left, 4: top left
input_parameters.backgroundEdge   = 100;
input_parameters.backgroundDist   = 0;

% configuration of correlation-based drift correction
input_parameters.correlationFlag  = 1;       % flag for enabling/disabling correlation-based drift correction

% configuration of global drift correction
input_parameters.globalMode       = 1;       % 0: disable global drift correction
                                             % 1: automatic global drift correction based on geometrical center computation
                                             % 2: manual global correction using vectors provided for reference time points

% parameters required when globalMode == 1
input_parameters.maskFactor       = 0.4;
input_parameters.scaling          = 2.031 / (6.5 / 16);
input_parameters.smoothing        = [1 20];  % smoothing flag ('rloess'), smoothing window size
input_parameters.kernelSize       = 51;
input_parameters.kernelSigma      = 20;
input_parameters.maskMinimum      = 1;       % percentile for minimum calculation (0 = true minimum)
input_parameters.fraction         = 0.01;    % minimal object size in binary mask, set to 0 to disable bwareaopen (default value 10 ^ -5)

% parameters required when globalMode == 2
input_parameters.referenceDrift   = [...     % each row follows this structure: reference time point (slot 1), x-/y-/z-drift vector (slots 2-4), follows Matlab convention
      0,  0,  0,  0; ...    % Note: at least two rows are required to enable manual global drift correction
    100, 10, 20, 30];

input_parameters.maxStampDigits   = 6;
input_parameters.poolWorkers      = 0;       % use "0" to enable automated detection of available CPU cores
%%
fn = [pwd '/sampleInputLocalEC.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));