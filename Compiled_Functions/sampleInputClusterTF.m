%% parameters
input_parameters = [];

input_parameters.timepoints.start = 0;
input_parameters.timepoints.every = 1;
input_parameters.timepoints.end = 50;

input_parameters.references.start   = 0;
input_parameters.references.every = 1;
input_parameters.references.end = 50;

input_parameters.inputString  = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected';
input_parameters.sourceString = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.MultiFused/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039';
input_parameters.outputString = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.TimeFused/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039';
input_parameters.inputID      = '_blending';
input_parameters.outputID     = '_blending';
input_parameters.lookUpTable  = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Configurations_Mmu/SPM00_CM00_CM01_CHN00_CHN01_analyzeParameters/lookUpTable.mat';
input_parameters.dataType     = 1;         % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks

input_parameters.specimen     = 0;
input_parameters.cameras      = 0:1;
input_parameters.sChannels    = 0:1;       % source channels
input_parameters.tChannels    = 0:1;       % target channels

input_parameters.reducedIO    = 1;         % 0 for full logging, 1 for minimal logging, 2 for slice logging
input_parameters.inputType    = 0;         % 0: input data in KLB format
                          % 1: input data in JP2 format
                          % 2: input data in TIF format
input_parameters.outputType   = 0;         % 0: output data saved in KLB format
                          % 1: output data saved in JP2 format
                          % 2: output data saved in TIF format
input_parameters.splitting    = 10;
input_parameters.intSizes     = [10 5];    % adaptive slab size for intensity correction, set to 0 to use slice mode (first entry: channels, second entry: cameras)
input_parameters.correction   = [3 3 0 0]; % slots 1 and 2 provide intensity correction flags for multi-channel and multi-camera: 0 to disable, 1 for calculated values, 2 for manual values (slots 3 and 4), 3 for lookUpTable
                          % slots 3 and 4 provide manual correction factors for the transformed channel/camera with respect to the reference channel/camera (will be inversely applied if < 1)
input_parameters.percentile   = 5;         % percentile for determining background levels
input_parameters.subSampling  = 100;       % data set subsampling for percentile computation
input_parameters.fusionType   = 0;         % 0 for blending, 1 for wavelet fusion, 2 for averaging
input_parameters.blending     = [20 4];    % blending ranges (first entry: channels, second entry: cameras)

input_parameters.enforceFlag  = [1 1, 1 1 1, 1, 1 1 1];
input_parameters.verbose      = 0;         % 0 for minimal text output, 1 for timing text output, 2 for full text ouput

input_parameters.cropping     = {[0 0 0 0 0 0]; [0 0 0 0 0 0]}; % ImageJ xy-ROI, start plane, stop plane (camera 0/1); set all to 0 to disable cropping
input_parameters.scaling      = 2.031 / (6.5 / 16);             % axial step size <divided by> (pixel pitch <divided by> magnification)

input_parameters.leftFlags    = [2 1];     % indicates in which channel the light sheet comes from the left-hand side (camera 0/1), 1 = reference channel, 2 = transformed channel

input_parameters.flipHFlag    = 1;         % indicates whether the stack recorded with the second camera should be flipped horizontally
input_parameters.flipVFlag    = 0;         % indicates whether the stack recorded with the second camera should be flipped vertically

input_parameters.frontFlag    = 1;         % indicates for which camera the high-quality information is in the front, 1 = reference camera, 2 = transformed camera

input_parameters.localRun     = [1 0];     % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                          %         note: cluster submission requires a Windows cluster with HPC 2008 support (see "job submit" commands below)
                          % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                          %         note: use "0" to enable automated detection of available CPU cores

input_parameters.jobMemory    = [1 0];     % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
                          % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                          %         note 1: slot 2 is only evaluated if automated memory management is disabled
                          %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

input_parameters.coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)
                                                       % note: parameter is only required for cluster submission
                          % DGA added globalMask                             
input_parameters.globalMask   = [0 0];     % slot 1 provides flag for use of a global (constant) mask, slot 2 provides reference time point for global mask
%%
fn = [pwd '/sampleInput_clusterTF.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));                          