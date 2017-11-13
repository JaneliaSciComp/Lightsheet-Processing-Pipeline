%% parameter settings

timepoints = [100 120 180 240 250 300 360 400 420 480]; % [60 100 120 180 240 250 300 360 400 420 480];

imageFileNameSegment1 = ['X:' filesep 'SV1' filesep '14-05-21' filesep 'Mmu_E1_CAGTAG1.corrected' filesep 'Results' filesep 'TimeFused.Corrected' filesep 'Mmu_E1_CAGTAG1.TM'];
imageFileNameSegment2 = ['_timeFused_blending' filesep 'SPM00_TM'];
imageFileNameSegment3 = '_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';

predictionFileNameSegment1  = ['X:' filesep 'SV1' filesep '14-05-21' filesep 'DivisionDetection' filesep 'Kristin' filesep 'MK5 Sparse' filesep];
predictionFileNameSegment2  = '.h5';
predictionsDataName         = '/coo';
writePredictionsSingleFlag  = 0;
constrainXYToPredictionData = 1;
constrainZToPredictionData  = 1;

annotationFileName    = 'annotations.mat';
annotationLateralSize = 10;
annotationAxialSize   = 5;

outputFolder = ['X:' filesep 'SV1' filesep '14-05-21' filesep 'DivisionDetection' filesep 'Kristin' filesep 'MK5 Visualization'];

for t = 1:numel(timepoints)
    timepoint = timepoints(t);
    
    disp(['processing timepoint ' num2str(timepoint)]);
    
    %% process image data
    
    if exist(outputFolder, 'dir') ~= 7
        mkdir(outputFolder);
    end;
    
    imageFileName = [imageFileNameSegment1 num2str(timepoint, '%.6d') imageFileNameSegment2 num2str(timepoint, '%.6d') imageFileNameSegment3];
    imageStack = permute(readImage(imageFileName), [2 1 3]);
    
    [xSize, ySize, zSize] = size(imageStack);
    
    %% process prediction data
    
    predictionsFileName = [predictionFileNameSegment1 num2str(timepoint) predictionFileNameSegment2];
    predictions = h5read(predictionsFileName, predictionsDataName);
    predictionCoordinates = double(predictions([2 3 1], :));
    predictionScores = predictions(4, :);
    boundaries = [...
        min(predictionCoordinates(1, :)), max(predictionCoordinates(1, :)), ...
        min(predictionCoordinates(2, :)), max(predictionCoordinates(2, :)), ...
        min(predictionCoordinates(3, :)), max(predictionCoordinates(3, :))];
    clear predictions;
    
    if constrainXYToPredictionData
        xStart = boundaries(1);
        xStop = boundaries(2);
        yStart = boundaries(3);
        yStop = boundaries(4);
    else
        xStart = 1;
        xStop = xSize;
        yStart = 1;
        yStop = ySize;
    end;
    
    if constrainZToPredictionData
        zStart = boundaries(5);
        zStop = boundaries(6);
    else
        zStart = 1;
        zStop = zSize;
    end;
    
    writeImage(imageStack(xStart:xStop, yStart:yStop, zStart:zStop), [outputFolder filesep 'Timepoint_' num2str(timepoint) '.Images.klb']);
    
    predictionSlots = sub2ind([xSize, ySize, zSize], predictionCoordinates(1, :), predictionCoordinates(2, :), predictionCoordinates(3, :));
    predictionStack = zeros([xSize, ySize, zSize], 'single');
    predictionStack(predictionSlots) = predictionScores;
    if writePredictionsSingleFlag
        writeImage(predictionStack(xStart:xStop, yStart:yStop, zStart:zStop), [outputFolder filesep 'Timepoint_' num2str(timepoint) '.Predictions.SinglePrecision.klb']);
    end;
    
    predictionStack = uint16(predictionStack .* (2^16 - 1));
    writeImage(predictionStack(xStart:xStop, yStart:yStop, zStart:zStop), [outputFolder filesep 'Timepoint_' num2str(timepoint) '.Predictions.klb']);
    
    %% process annotation data
    
    load(annotationFileName);
    annotations = annotations(annotations(:, 1) == timepoint, 2:4);
    annotations(:, 3) = round(annotations(:, 3) ./ 5);
    annotations = uint16(round(annotations(:, [2 1 3])));
    
    annotationStack = zeros([xSize, ySize, zSize], 'uint16');
    for i = 1:size(annotations, 1)
        currentLateralMarkerSize = min([annotationLateralSize, ...
            annotations(i, 1) - 1, size(imageStack, 1) - annotations(i, 1), ...
            annotations(2, 1) - 1, size(imageStack, 2) - annotations(i, 2)]);
        y = -double(currentLateralMarkerSize):0.1:double(currentLateralMarkerSize);
        yCoordinates = double(annotations(i, 2)) + round(y);
        xCoordinatesPlus = double(annotations(i, 1)) + round(sqrt(double(currentLateralMarkerSize) ^ 2 - y .^ 2));
        xCoordinatesMinus = double(annotations(i, 1)) - round(sqrt(double(currentLateralMarkerSize) ^ 2 - y .^ 2));
        for z = max(1, annotations(i, 3) - annotationAxialSize):min(size(imageStack, 3), annotations(i, 3) + annotationAxialSize)
            for j = 1:numel(yCoordinates)
                annotationStack(xCoordinatesPlus(j), yCoordinates(j), z) = 2^16 - 1;
                annotationStack(xCoordinatesMinus(j), yCoordinates(j), z) = 2^16 - 1;
            end;
        end;
        
        currentLateralMarkerSize = currentLateralMarkerSize - 1;
        y = -double(currentLateralMarkerSize):0.1:double(currentLateralMarkerSize);
        yCoordinates = double(annotations(i, 2)) + round(y);
        xCoordinatesPlus = double(annotations(i, 1)) + round(sqrt(double(currentLateralMarkerSize) ^ 2 - y .^ 2));
        xCoordinatesMinus = double(annotations(i, 1)) - round(sqrt(double(currentLateralMarkerSize) ^ 2 - y .^ 2));
        for z = max(1, annotations(i, 3) - annotationAxialSize):min(size(imageStack, 3), annotations(i, 3) + annotationAxialSize)
            for j = 1:numel(yCoordinates)
                annotationStack(xCoordinatesPlus(j), yCoordinates(j), z) = 2^16 - 1;
                annotationStack(xCoordinatesMinus(j), yCoordinates(j), z) = 2^16 - 1;
            end;
        end;
    end;
    writeImage(annotationStack(xStart:xStop, yStart:yStop, zStart:zStop), [outputFolder filesep 'Timepoint_' num2str(timepoint) '.Annotations.klb']);
end;