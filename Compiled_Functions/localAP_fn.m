function localAP_fn(filename)
input_parameters = loadjson(fileread(filename));
input_parameters = convert_limits_to_values(input_parameters, {'timepoints','fullInterval'});
if ~isfield(input_parameters, 'verbose'), input_parameters.verbose = true; end
if input_parameters.verbose
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end
timepoints   = input_parameters.timepoints;
fullInterval = input_parameters.fullInterval;

outputString = input_parameters.outputString;
outputID     = input_parameters.outputID;
dataType     = input_parameters.dataType;      % 0 for unsegmented clusterPT output stacks, 1 for segmented clusterPT output stacks

specimen     = input_parameters.specimen;
angle        = input_parameters.angle;
cameras      = input_parameters.cameras;
channels     = input_parameters.channels;

readFactors  = input_parameters.readFactors;      % 0 to ignore intensity correction data, 1 to include intensity correction data
smoothing    = input_parameters.smoothing; % smoothing flag ('rloess'), smoothing window size
offsetRange  = input_parameters.offsetRange;     % averaging range for offsets
angleRange   = input_parameters.angleRange;     % averaging range for angles
intRange     = input_parameters.intRange;      % averaging range for intensity correction
averaging    = input_parameters.averaging;      % 0 for mean, 1 for median
staticFlag   = input_parameters.staticFlag;      % 0 for time-dependent registration and intensity correction, 1 for static parameter set

configRoot   = input_parameters.configRoot;

analyzeParameters(...
    timepoints, fullInterval,...
    outputString, outputID, dataType,...
    specimen, angle, cameras, channels, readFactors, smoothing,...
    offsetRange, angleRange, intRange, averaging, staticFlag, configRoot)
end