input_parameters.timepoints   = 0:50;
input_parameters.fullInterval = 0:50;

input_parameters.outputString = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/from_kellerp/Example Data/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039.correctedMultiFused/Mmu_E1_H2BeGFPxTmCherry_01_20170928_144039';
input_parameters.outputID     = '_blending';
input_parameters.dataType     = 1;      % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks

input_parameters.specimen     = 0;
input_parameters.cameras      = 0:1;
input_parameters.channels     = 0:1;

input_parameters.readFactors  = 1;      % 0 to ignore intensity correction data, 1 to include intensity correction data
input_parameters.smoothing    = [1 20]; % smoothing flag ('rloess'), smoothing window size
input_parameters.offsetRange  = 10;     % averaging range for offsets
input_parameters.angleRange   = 10;     % averaging range for angles
input_parameters.intRange     = 5;      % averaging range for intensity correction
input_parameters.averaging    = 0;      % 0 for mean, 1 for median
                       % DGA added
input_parameters.staticFlag   = 0;      % 0 for time-dependent registration and intensity correction, 1 for static parameter set
%%
fn = [pwd '/sampleLocalAPInput.json'];
str = savejson('', input_parameters);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
parameters = loadjson(fileread(fn));