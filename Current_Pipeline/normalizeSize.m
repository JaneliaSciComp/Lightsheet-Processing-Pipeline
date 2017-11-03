%% parameters

inputRoot       = 'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Results' filesep 'MultiFused';
outputRoot      = 'X:' filesep 'SiMView1' filesep '14-01-21' filesep 'Mmu_E1_CAGTAG1_01_23_20140121_141339.corrected' filesep 'Results' filesep 'MultiFused.NormalizedSize';
headerPattern   = 'Mmu_E1_CAGTAG1.TM??????_multiFused_blending' filesep;
filePattern     = 'SPM00_TM??????_CM00_CM01_CHN00_CHN01.fusedStack';

timepoints      = 0:305;
percentile      = [1 10]; % slot 1: percentile for minimum calculation (0 = true minimum); slot 2: subsampling
target          = 2;      % 0: stacks only, 1: projections only, 2: stacks and projections

inputType       = 0;      % 0: input data in KLB format
                          % 1: input data in JP2 format
                          % 2: input data in TIF format

maxStampDigits  = 6;
poolWorkers     = 0;      % use "0" to enable automated detection of available CPU cores

%% main loop

if poolWorkers == 0
    poolWorkers = feature('numcores');
    disp(' ');
    disp([num2str(poolWorkers) ' CPU cores were detected and will be allocated for parallel processing.']);
end;

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end;

nTimepoints = numel(timepoints);
inputDatabase = cell(nTimepoints, 4);
outputDatabase = cell(nTimepoints, 4);

for n = 1:nTimepoints
    outputDatabase{n, 1} = [outputRoot filesep '' headerPattern filePattern];
    for i = maxStampDigits:-1:1
        positions = strfind(outputDatabase{n, 1}, repmat('?', [1 i]));
        precision = ['%.' num2str(i) 'd'];
        for k = 1:length(positions)
            outputDatabase{n, 1}(positions(k):(positions(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
    end;
    
    outputDatabase{n, 5} = [outputRoot filesep '' headerPattern];
    for i = maxStampDigits:-1:1
        positions = strfind(outputDatabase{n, 5}, repmat('?', [1 i]));
        precision = ['%.' num2str(i) 'd'];
        for k = 1:length(positions)
            outputDatabase{n, 5}(positions(k):(positions(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
    end;
    
    outputDatabase{n, 2} = [outputDatabase{n, 1} '_xyProjection' inputExtension];
    outputDatabase{n, 3} = [outputDatabase{n, 1} '_xzProjection' inputExtension];
    outputDatabase{n, 4} = [outputDatabase{n, 1} '_yzProjection' inputExtension];
    
    outputDatabase{n, 1} = [outputDatabase{n, 1} inputExtension];
    
    inputDatabase{n, 1} = [inputRoot filesep '' headerPattern filePattern];
    for i = maxStampDigits:-1:1
        positions = strfind(inputDatabase{n, 1}, repmat('?', [1 i]));
        precision = ['%.' num2str(i) 'd'];
        for k = 1:length(positions)
            inputDatabase{n, 1}(positions(k):(positions(k) + i - 1)) = num2str(timepoints(n), precision);
        end;
    end;
    
    inputDatabase{n, 2} = [inputDatabase{n, 1} '_xyProjection' inputExtension];
    inputDatabase{n, 3} = [inputDatabase{n, 1} '_xzProjection' inputExtension];
    inputDatabase{n, 4} = [inputDatabase{n, 1} '_yzProjection' inputExtension];
    
    inputDatabase{n, 1} = [inputDatabase{n, 1} inputExtension];
end;

disp(' ');
disp('Collecting information on image volume size as a function of time');

dimensions = zeros(numel(timepoints), 4);
dimensions(:, 1) = timepoints;

for n = 1:nTimepoints
    switch inputType
        case 0
            headerInformation = readKLBheader(inputDatabase{n, 1});
            stackDimensions = headerInformation.xyzct(1:3);
        case 1
            [stackDimensions, bitDepth] = readJP2header(inputDatabase{n, 1});
        case 2
            headerInformation = imfinfo(inputDatabase{n, 1});
            stackDimensions = [headerInformation(1).Height headerInformation(1).Width numel(headerInformation)];
    end;
    
    dimensions(n, 2) = stackDimensions(1);
    dimensions(n, 3) = stackDimensions(2);
    dimensions(n, 4) = stackDimensions(3);
end;

save([outputRoot '\dimensions.mat'], 'dimensions');

maxDimensions = max(dimensions(:, 2:4), [], 1);
save([outputRoot filesep 'maxDimensions.mat'], 'maxDimensions');

deltas = zeros(numel(timepoints), 4);
deltas(:, 1) = timepoints;

for n = 1:nTimepoints
    deltas(n, 2) = floor((maxDimensions(1) - dimensions(n, 2)) / 2);
    deltas(n, 3) = floor((maxDimensions(2) - dimensions(n, 3)) / 2);
    deltas(n, 4) = floor((maxDimensions(3) - dimensions(n, 4)) / 2);
end;

save([outputRoot '\deltas.mat'], 'deltas');

disp(' ');
disp(['Parallelizing size normalization for ' num2str(nTimepoints) ' time points.']);
disp(' ');

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(poolWorkers);
disp(' ');

parfor n = 1:nTimepoints
    switch target
        case 0
            disp(['Normalizing size of image stack at time point ' num2str(timepoints(n), '%.6d')]);
        case 1
            disp(['Normalizing size of stack projections at time point ' num2str(timepoints(n), '%.6d')]);
        case 2
            disp(['Normalizing size of image stack and stack projections at time point ' num2str(timepoints(n), '%.6d')]);
    end;
    
    if ~exist(outputDatabase{n, 5}, 'dir')
        mkdir(outputDatabase{n, 5});
    end;
    
    dimensionsIndex = find(dimensions(:, 1) == timepoints(n), 1);
    mismatchFlag = ...
        (maxDimensions(1) ~= dimensions(dimensionsIndex, 2)) || ...
        (maxDimensions(2) ~= dimensions(dimensionsIndex, 3)) || ...
        (maxDimensions(3) ~= dimensions(dimensionsIndex, 4));
    
    if target == 0 || target == 2
        if mismatchFlag
            unpaddedStack = readImage(inputDatabase{n, 1});
            
            if percentile(1) == 0
                currentStackMin = min(unpaddedStack(1:percentile(2):end));
            else
                currentStackMin = prctile(unpaddedStack(1:percentile(2):end), percentile(1));
            end;
            currentStack = ones(maxDimensions(1), maxDimensions(2), maxDimensions(3), 'uint16') .* currentStackMin;
            currentStack(...
                (deltas(dimensionsIndex, 2) + 1):(deltas(dimensionsIndex, 2) + dimensions(dimensionsIndex, 2)), ...
                (deltas(dimensionsIndex, 3) + 1):(deltas(dimensionsIndex, 3) + dimensions(dimensionsIndex, 3)), ...
                (deltas(dimensionsIndex, 4) + 1):(deltas(dimensionsIndex, 4) + dimensions(dimensionsIndex, 4))) = unpaddedStack;
            
            writeImage(currentStack, outputDatabase{n, 1});
        else
            copyfile(inputDatabase{n, 1}, outputDatabase{n, 1});
        end;
    end;
    
    if target == 1 || target == 2
        if mismatchFlag
            unpaddedImageXY = readImage(inputDatabase{n, 2});
            unpaddedImageXZ = readImage(inputDatabase{n, 3});
            unpaddedImageYZ = readImage(inputDatabase{n, 4});
            
            if percentile(1) == 0
                currentImageXYMin = min(unpaddedImageXY(:));
                currentImageXZMin = min(unpaddedImageXZ(:));
                currentImageYZMin = min(unpaddedImageYZ(:));
            else
                currentImageXYMin = prctile(unpaddedImageXY(:), percentile(1));
                currentImageXZMin = prctile(unpaddedImageXZ(:), percentile(1));
                currentImageYZMin = prctile(unpaddedImageYZ(:), percentile(1));
            end;
            
            currentImageXY = ones(maxDimensions(1), maxDimensions(2), 'uint16') .* currentImageXYMin;
            currentImageXZ = ones(maxDimensions(1), maxDimensions(3), 'uint16') .* currentImageXZMin;
            currentImageYZ = ones(maxDimensions(2), maxDimensions(3), 'uint16') .* currentImageYZMin;
            
            currentImageXY((deltas(n, 2) + 1):(deltas(n, 2) + dimensions(n, 2)), (deltas(n, 3) + 1):(deltas(n, 3) + dimensions(n, 3))) = unpaddedImageXY;
            currentImageXZ((deltas(n, 2) + 1):(deltas(n, 2) + dimensions(n, 2)), (deltas(n, 4) + 1):(deltas(n, 4) + dimensions(n, 4))) = unpaddedImageXZ;
            currentImageYZ((deltas(n, 3) + 1):(deltas(n, 3) + dimensions(n, 3)), (deltas(n, 4) + 1):(deltas(n, 4) + dimensions(n, 4))) = unpaddedImageYZ;
            
            writeImage(currentImageXY, outputDatabase{n, 2});
            writeImage(currentImageXZ, outputDatabase{n, 3});
            writeImage(currentImageYZ, outputDatabase{n, 4});
        else
            copyfile(inputDatabase{n, 2}, outputDatabase{n, 2});
            copyfile(inputDatabase{n, 3}, outputDatabase{n, 3});
            copyfile(inputDatabase{n, 4}, outputDatabase{n, 4});
        end;
    end;
end;

disp(' ');

if matlabpool('size') > 0
    matlabpool('close');
end;

disp(' ');