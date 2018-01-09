%% keller parameters
input_parameters = [];
input_parameters.inputRoot        = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.TimeFused/';
input_parameters.outputRoot       = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.CorrectStack/';
input_parameters.headerPattern    = 'Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.TM??????_timeFused_blending/';
input_parameters.filePattern      = 'SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';
input_parameters.configRoot       = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Configurations/SPM00_CM00_CM01_CHN00_CHN01_stackCorrection';

input_parameters.timepoints       = 0:50;
input_parameters.dataType         = 1;      % 0: process projections only, 1: process stacks and projections
input_parameters.percentile       = [1 10]; % [percentileFlag, subsampling], percentileFlag - 0: calculate true minimum, otherwise: use as percentile for minimum calculation

input_parameters.inputType        = 0;      % 0: input data in KLB format
                                            % 1: input data in JP2 format
                                            % 2: input data in TIF format
input_parameters.outputType       = 0;      % 0: output data saved in KLB format
                                            % 1: output data saved in JP2 format
                                            % 2: output data saved in TIF format

% configuration of drift correction
input_parameters.correctDrift     = 1;
input_parameters.referenceTime    = 25;
input_parameters.referenceROI     = [];     % [xStart xStop; yStart yStop; zStart zStop], provide empty vector to enforce use of maximum dimensions

% configuration of intensity normalization
input_parameters.correctIntensity = 1;

input_parameters.maxStampDigits   = 6;

input_parameters.localRun         = [0 0];  % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                                            % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                                            %         note: use "0" to enable automated detection of available CPU cores

input_parameters.jobMemory        = [1 0];  % slot 1: flag for automated memory management (0: disable, 1: enable)
                                            % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                                            %         note 1: slot 2 is only evaluated if automated memory management is disabled
                                            %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

input_parameters.coreMemory       = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)
input_parameters.verbose      = 0; 
%%
fn = [pwd '/sampleInput_clusterCS.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));