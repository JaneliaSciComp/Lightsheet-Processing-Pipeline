%% configuration

logFilename  = ['T:' filesep 'SiMView1' filesep '15-04-21' filesep 'Scripts' filesep 'createMasks.processingLogs_SPM0.mat'];
              %'T:' filesep 'SiMView1' filesep '15-04-21' filesep 'Scripts' filesep 'createMasks.processingLogs_SPM1.mat';

dataFolder   = ['S:' filesep 'SiMView1' filesep '15-04-21' filesep 'Mmu_E1_mKate2_20150421_143006.corrected'];
specimen     = 0;
              %1;

timepoints   = 0:1243;
channel      = 0;
cameras      = 0:1;

identifier   = 'unmasked_r5_d3_e2.threshold_0.2'; % identification string used in file names for calculateWeights output data

outputFolder = ['T:' filesep 'SiMView1' filesep '15-04-21' filesep 'Scripts' filesep 'WeightedStackProjections.SPM0'];
              %'T:' filesep 'SiMView1' filesep '15-04-21' filesep 'Scripts' filesep 'WeightedStackProjections.SPM1';

background   = [1 100];

%% main loop

if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end;

if ~isempty(identifier)
    identifier = ['_' identifier];
end;

if ~isempty(logFilename)
    load(logFilename, 'processedMetaData'); % timepoints (1), zLowPlaneAdditions (5), zHighPlaneAdditions (6)
    
    zLowPlaneAdditions = [];
    zHighPlaneAdditions = [];
    for t = 1:numel(timepoints)
        zLowPlaneAdditions = cat(1, zLowPlaneAdditions, processedMetaData(processedMetaData(:, 1) == timepoints(t), 5));
        zHighPlaneAdditions = cat(1, zHighPlaneAdditions, processedMetaData(processedMetaData(:, 1) == timepoints(t), 6));
    end;
end;

xyHeader = readKLBheader([dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(end), '%.6d') ...
    filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(end), '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' identifier '.xyProjection.klb']);
xzHeader = readKLBheader([dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(end), '%.6d') ...
    filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(end), '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' identifier '.xzProjection.klb']);
xSize = xyHeader.xyzct(1);
ySize = xyHeader.xyzct(2);
zSize = xzHeader.xyzct(2);

for camera = cameras
    disp(' ');
    disp(['Processing camera ' num2str(camera)]);
    
    if background(1)
        xySequence = ones(xSize, ySize, numel(timepoints), 'uint16').*background(2);
        xzSequence = ones(xSize, zSize, numel(timepoints), 'uint16').*background(2);
        yzSequence = ones(ySize, zSize, numel(timepoints), 'uint16').*background(2);
    else
        xySequence = zeros(xSize, ySize, numel(timepoints), 'uint16');
        xzSequence = zeros(xSize, zSize, numel(timepoints), 'uint16');
        yzSequence = zeros(ySize, zSize, numel(timepoints), 'uint16');
    end;    
    
    for t = 1:numel(timepoints)
        disp(['* Integrating time point ' num2str(timepoints(t))]);
        
        weightedStackXYProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' identifier '.xyProjection.klb'];
        weightedStackXZProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' identifier '.xzProjection.klb'];
        weightedStackYZProjectionFilename = [dataFolder filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoints(t), '%.6d') ...
            filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoints(t), '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channel, '%.2d') '.weighted' identifier '.yzProjection.klb'];
        xyProjection = readImage(weightedStackXYProjectionFilename);
        xzProjection = readImage(weightedStackXZProjectionFilename);
        yzProjection = readImage(weightedStackYZProjectionFilename);
        
        xySequence(:, :, t) = xyProjection;
        if ~isempty(logFilename)
            xzSequence(:, (1+zLowPlaneAdditions(end)-zLowPlaneAdditions(t)):(1+zLowPlaneAdditions(end)-zLowPlaneAdditions(t)+size(xzProjection, 2)-1), t) = xzProjection;
            yzSequence(:, (1+zLowPlaneAdditions(end)-zLowPlaneAdditions(t)):(1+zLowPlaneAdditions(end)-zLowPlaneAdditions(t)+size(yzProjection, 2)-1), t) = yzProjection;
        else
            xzSequence(:, (1+floor((zSize-size(xzProjection, 2))/2)):(1+floor((zSize-size(xzProjection, 2))/2)+size(xzProjection, 2)-1), t) = xzProjection;
            yzSequence(:, (1+floor((zSize-size(yzProjection, 2))/2)):(1+floor((zSize-size(yzProjection, 2))/2)+size(yzProjection, 2)-1), t) = yzProjection;
        end;
    end;
    
    disp(' ');
    disp('Saving assembled image sequences');
    
    assembledXYSequenceFilename = [outputFolder filesep 'assembledWeightedProjections_CM' num2str(camera) '.xy.klb'];
    assembledXZSequenceFilename = [outputFolder filesep 'assembledWeightedProjections_CM' num2str(camera) '.xz.klb'];
    assembledYZSequenceFilename = [outputFolder filesep 'assembledWeightedProjections_CM' num2str(camera) '.yz.klb'];
    writeImage(xySequence, assembledXYSequenceFilename);
    writeImage(xzSequence, assembledXZSequenceFilename);
    writeImage(yzSequence, assembledYZSequenceFilename);
end;

disp(' ');