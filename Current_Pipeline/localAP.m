timepoints   = 0:0;
fullInterval = 0:0;

outputString = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/Paper/Supplementary_Data_2/Image_Data.MultiFused/Dme_E1_H2ARFP';
outputID     = '_blending';
configRoot = '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet_Data/Paper/Supplementary_Data_2/configRoot';
dataType     = 1;      % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks

specimen     = 0;
cameras      = 0:1;
channels     = 0:1;

readFactors  = 1;      % 0 to ignore intensity correction data, 1 to include intensity correction data
smoothing    = [1 20]; % smoothing flag ('rloess'), smoothing window size
offsetRange  = 10;     % averaging range for offsets
angleRange   = 10;     % averaging range for angles
intRange     = 5;      % averaging range for intensity correction
averaging    = 0;      % 0 for mean, 1 for median
                       % DGA added
staticFlag   = 0;      % 0 for time-dependent registration and intensity correction, 1 for static parameter set

analyzeParameters(...
    timepoints, fullInterval,...
    outputString, outputID, dataType,...
    specimen, cameras, channels, readFactors, smoothing,...
    offsetRange, angleRange, intRange, averaging, staticFlag, configRoot)
%% from keller
% % % timepoints   = 0:5:733;
% % % fullInterval = 0:733;
% % % 
% % % outputString = 'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'MultiFused' filesep 'Pha_E1_H2bRFP';
% % % outputID     = '_blending';
% % % dataType     = 1;      % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks
% % % 
% % % specimen     = 0;
% % % cameras      = 0:1;
% % % channels     = 2:3;
% % % 
% % % readFactors  = 1;      % 0 to ignore intensity correction data, 1 to include intensity correction data
% % % smoothing    = [1 20]; % smoothing flag ('rloess'), smoothing window size
% % % offsetRange  = 10;     % averaging range for offsets
% % % angleRange   = 10;     % averaging range for angles
% % % intRange     = 5;      % averaging range for intensity correction
% % % averaging    = 0;      % 0 for mean, 1 for median
% % % staticFlag   = 0;      % 0 for time-dependent registration and intensity correction, 1 for static parameter set
% % % 
% % % analyzeParameters(...
% % %     timepoints, fullInterval,...
% % %     outputString, outputID, dataType,...
% % %     specimen, cameras, channels, readFactors, smoothing,...
% % %     offsetRange, angleRange, intRange, averaging, staticFlag)