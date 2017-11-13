clear all;
input_parameters.inputFolder  = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/Paper/Supplementary_Data_1/Image_Data';
input_parameters.outputLabel  = '';
input_parameters.specimen     = 0;                   % specimen index to be processed
input_parameters.timepoints   = 0;                   % time points to be processed
input_parameters.cameras      = [0 1];               % camera indices to be processed
input_parameters.channels     = [0 1];               % channel indices to be processed

input_parameters.dimensions   = [603 1272 125];      % override parameter for user-provided raw stack dimensions (ImageJ-x, ImageJ-y, z/t)
                                    % note: assign empty vector [] to indicate that override is not needed (requires XML meta data)

% use highly saturated global XY projection to determine the following four parameters (can be left empty to enforce maximum ROI)
input_parameters.startsLeft   = [  75   24];         % cropping left start coordinates for each camera (ImageJ-x convention)
input_parameters.startsTop    = [  24   17];         % cropping top start coordinates for each camera (ImageJ-y convention)
input_parameters.widths       = [ 500  500];         % cropping width for each camera (ImageJ-w convention)
input_parameters.heights      = [1172 1172];         % cropping height for each camera (ImageJ-h convention)

% use highly saturated global XZ projection to determine the following two parameters (can be left empty to enforce maximum ROI)
input_parameters.startsFront  = [   3    3];         % cropping front start coordinates for each camera (ImageJ-y convention)
input_parameters.depths       = [ 120  120];         % cropping depth for each camera (ImageJ-h convention)

input_parameters.inputType    = 0;                   % 0: input data in TIF format
                                    % 1: input data in JP2 format
                                    % 2: input data in binary stack format (normal experiment mode, structured)
                                    % 3: input data in binary stack format (normal experiment mode, unstructured)
                                    % 4: input data in binary stack format (HS single-plane experiment mode)
                                    %    note: pixel correction and segmentFlag are inactive in this mode
input_parameters.outputType   = 0;                   % 0: output data saved in KLB format
                                    % 1: output data saved in JP2 format
                                    % 2: output data saved in TIF format
input_parameters.correctTIFF  = 0;                   % 0: never transpose output data, 1: transpose output data if TIFF parity is odd

input_parameters.rotationFlag = 0;                   % 0: do not rotate image stacks, 1: rotate image stacks by 90 degrees clockwise, -1: rotate image stacks by 90 degrees counter-clockwise
input_parameters.medianRange  = [3 3];               % kernel x/y-size for median filter for dead pixel detection and removal
input_parameters.percentile   = [1 5 100];           % slot 1: background percentile for mask calculation
                                    % slot 2: background percentile for intensity correction
                                    % slot 3: volume sub-sampling for background estimation

input_parameters.segmentFlag  = 1;                   % 0: do not remove background in output stacks, 1: remove background in output stacks (default)
input_parameters.flipHFlag    = 0;                   % indicates whether the stack recorded with the second camera should be flipped horizontally
input_parameters.flipVFlag    = 0;                   % indicates whether the stack recorded with the second camera should be flipped vertically
input_parameters.splitting    = 10;                  % level of stack splitting when performing Gauss convolution
input_parameters.kernelSize   = 5;                   % Gauss kernel size
input_parameters.kernelSigma  = 2;                   % Gauss kernel sigma
input_parameters.scaling      = 2.031 / (6.5 / 16);  % axial step size <divided by> (pixel pitch <divided by> magnification)
input_parameters.references   = [0 1];               % list of reference channel groups that define segmentation masks for (content-)dependent channel groups
                                    % note 1: leave empty if all channels are to be treated as independent channels
                                    % note 2: provide multiple elements per row (i.e. a channel group) to fuse segmentation masks
input_parameters.dependents   = [];                  % list of (content-)dependent channel groups that relate to each reference channel group
input_parameters.thresholds   = [0.5 0.5];           % adaptive thresholds for each reference channel group
                                    % note: provide single threshold if "references" is left empty

input_parameters.loggingFlag  = 0;                   % 0: do not save processing logs (default), 1: save processing logs
input_parameters.verbose      = 0;                   % 0: do not display processing information, 1: display processing information

input_parameters.localRun     = [1 0];               % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                                    %         note: cluster submission requires a Windows cluster with HPC 2008 support (see "job submit" commands below)
                                    % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                                    %         note: use "0" to enable automated detection of available CPU cores

input_parameters.jobMemory    = [1 0];               % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
                                    %         note: setting jobMemory(1) to "2" is incompatible with inputType == 4 or numel(dimensions) ~= 0
                                    % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                                    %         note 1: slot 2 is only evaluated if automated memory management is disabled
                                    %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

input_parameters.coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)
                                                       % note: parameter is only required for cluster submission

fn = [pwd '/sampleClusterPTInput.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));


