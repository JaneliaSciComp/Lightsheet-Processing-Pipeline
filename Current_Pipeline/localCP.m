timepoints     = 0:733;
configurations = {...
    'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'MultiFused' filesep 'Pha_E1_H2bRFP',          '_multiFused_blending', '',             0, 0:1, 0:1, 'Projections.MultiFused'; ...
    'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'TimeFused' filesep 'Pha_E1_H2bRFP',           '_timeFused_blending',  '',             0, 0:1, 0:1, 'Projections.TimeFused'; ...
    'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'TimeFused.Filtered' filesep 'Pha_E1_H2bRFP',  '_timeFused_blending',  'filtered_100', 0, 0:1, 0:1, 'Projections.TimeFused.Filtered'; ...
    'X:' filesep 'SiMView2' filesep '13-12-30' filesep 'Pha_E1_H2bRFP_01_20131230_140802.corrected' filesep 'Results' filesep 'TimeFused.Corrected' filesep 'Pha_E1_H2bRFP', '_timeFused_blending',  'corrected',    0, 0:1, 0:1, 'Projections.TimeFused.Corrected'; ...
    };

inputType    = 0; % 0: input data in KLB format
                  % 1: input data in JP2 format
                  % 2: input data in TIF format
outputType   = 0; % 0: output data saved in KLB format
                  % 1: output data saved in JP2 format
                  % 2: output data saved in TIF format

for i = 1:size(configurations, 1)
    outputString = configurations{i, 1};
    outputID     = configurations{i, 2};
    filterID     = configurations{i, 3};
    specimen     = configurations{i, 4};
    cameras      = configurations{i, 5};
    channels     = configurations{i, 6}; % processed channels (multiFuse, timeFuse)
    folder       = configurations{i, 7};

    collectProjections(...
        timepoints, outputString, outputID, filterID, ...
        specimen, cameras, channels, folder, ...
        inputType, outputType);
end;