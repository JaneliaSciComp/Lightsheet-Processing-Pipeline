%% keller parameters

input_parameters.timepoints   = 0:50;

input_parameters.inputDir     = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.CorrectStack/';
input_parameters.outputDir    = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.corrected.Filtered/';
input_parameters.header       = 'Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039';
input_parameters.footer       = '_timeFused_blending';
input_parameters.stackLabel   = '.corrected';

input_parameters.specimen     = 0;
input_parameters.cameras      = 0:1;
input_parameters.channels     = 0:1;

input_parameters.removeDirt   = [0 10]; % slot 1: flag for removing all but the largest connected component determined by thresholding
                                        % slot 2: threshold for binarizing volume

input_parameters.filterMode   = 2;      % [] = no filtering, 0 = subtract median-filtered image, 1 = subtract mean-filtered image, 2 = subtract Gauss-filtered image
input_parameters.rangeArray   = [100, 200];
input_parameters.splitting    = 15;
input_parameters.scaling      = 2.031 / (6.5 / 16);

input_parameters.preMedian    = [0 3];  % slot 1: flag for application of 2D median filter across image stack before background subtraction
                                        % slot 2: size of median filter kernel (resulting in command medfilt2(A, [n n]), where n is the kernel size)
                                        %         note: the filtering defined by filterMode will always be applied to the raw stack, so this is a separate processing branch
input_parameters.postMedian   = [0 3];  % slot 1: flag for application of 2D median filter across image stack after background subtraction
                                        % slot 2: size of median filter kernel (resulting in command medfilt2(A, [n n]), where n is the kernel size)
                                        %         note: if median filtering is activated, the result is stored separately from the background subtracted data

input_parameters.inputType    = 0;      % 0: input data in KLB format
                                        % 1: input data in JP2 format
                                        % 2: input data in TIF format
input_parameters.outputType   = 0;      % 0: output data saved in KLB format
                                        % 1: output data saved in JP2 format
                                        % 2: output data saved in TIF format
input_parameters.subProject   = 1;      % flag indicates whether separate front and back projections are created
input_parameters.saveRawMax   = 1;      % flag indicates whether projections of raw images are stored as well
input_parameters.saveStacks   = 1;      % flag indicates whether filtered image stacks are written to disk

input_parameters.localRun     = [0 0]; % slot 1: flag for local vs. cluster execution (0: cluster submission, 1: local workstation)
                                        % slot 2: number of parallel workers for execution on local workstation (only needed if slot 1 is set to 1)
                                        %         note: use "0" to enable automated detection of available CPU cores

input_parameters.jobMemory    = [0 0];  % slot 1: flag for automated memory management (0: disable, 1: enable, 2: enable and time-dependent)
                                        % slot 2: estimated upper boundary for memory consumption per submitted time point (in GB)
                                        %         note 1: slot 2 is only evaluated if automated memory management is disabled
                                        %         note 2: "0" indicates memory consumption below "coreMemory" threshold and enables parametric submission mode

input_parameters.coreMemory   = floor(((96 - 8) * 1024) / (12 * 1024)); % memory boundary for switching from parametric to memory-managed submission (in GB)
%%
input_parameters.verbose      = 0; 
fn = [pwd '/sampleInputClusterFR.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));