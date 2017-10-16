% annotations: t, x, y, z (isotropic)
% predictions: x, y, z (anisotropic), t, score
% thresholds:  threshold, test precision, test recall, train precision, train recall

testTimepoints    = [120, 240, 360];
trainTimepoints   = [100, 180, 250, 300, 400, 420 480];

searchRadius      = [20, 20, 10, 2]; % x, y, z (anisotropic), t
croppingRadius    = [80, 80, 10, 5]; % [40, 40, 10, 5]; % x, y, z (anisotropic), t
anisotropy        = 5;

FNSearchThreshold  = 0.1;
processAnnotations = 1;
processPredictions = 1;

dataFolder        = 'E:\Mouse Development\Division Detection\MK5 Sparse';

imagePathSegment1 = 'X:\SV1\14-05-21\Mmu_E1_CAGTAG1.corrected\Results\TimeFused.Corrected\Mmu_E1_CAGTAG1.TM';
imagePathSegment2 = '_timeFused_blending\SPM00_TM';
imagePathSegment3 = '_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';

outputFolder      = 'E:\Mouse Development\Division Detection\MK5 Inspection';

saveVolumes       = 0;
saveProjections   = 1;

%% processing predictions

if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end;

if exist([outputFolder '\Tables'], 'dir') ~= 7
    mkdir([outputFolder '\Tables']);
end;

load annotations;
annotations(:, 4) = annotations(:, 4) ./ anisotropy;

load([dataFolder '\predictions.mat']);
load([dataFolder '\thresholds.mat']);

disp(' ');

for c = 1:2
    switch c
        case 1
            timepoints = testTimepoints;
            displayLabel = 'testing';
            outputLabel = 'Testing';
        case 2
            timepoints = trainTimepoints;
            displayLabel = 'training';
            outputLabel = 'Training';
    end;
    
    annotatedPredictions = [];
    
    for t = timepoints
        annotatedPredictions = cat(1, annotatedPredictions, predictions(predictions(:, 4) == t, :));
    end;
    
    annotatedPredictions = cat(2, annotatedPredictions, zeros(size(annotatedPredictions, 1), 1));
    
    for t = timepoints
        extractionTimepoints = (t - croppingRadius(4)):(t + croppingRadius(4));
        for e = 1:numel(extractionTimepoints)
            disp(['reading image data at timepoint ' num2str(extractionTimepoints(e))]);
            
            stackPath = [imagePathSegment1 num2str(extractionTimepoints(e), '%.6d') ...
                imagePathSegment2 num2str(extractionTimepoints(e), '%.6d') imagePathSegment3];
            if e == 1
                stack = readImage(stackPath);
                imageBlock = zeros(size(stack, 1), size(stack, 2), size(stack, 3), numel(extractionTimepoints), 'uint16');
                imageBlock(:, :, :, 1) = stack;
            else
                imageBlock(:, :, :, e) = readImage(stackPath);
            end;
            
        end;
        
        if processAnnotations
            disp(['processing annotations for ' displayLabel ' timepoint ' num2str(t)]);
            
            annotationSelection = round(annotations(annotations(:, 1) == t, :));
            annotationSelection = cat(2, annotationSelection, zeros(size(annotationSelection, 1), 1));
            
            filteredPredictions = predictions(...
                predictions(:, 4) >= (t - searchRadius(4)) & predictions(:, 4) <= (t + searchRadius(4)) & ...
                predictions(:, 5) >= FNSearchThreshold, :);
            
            for i = 1:size(annotationSelection, 1)
                predictionDistances = filteredPredictions;
                predictionDistances(:, 1) = predictionDistances(:, 1) - annotationSelection(i, 2);
                predictionDistances(:, 2) = predictionDistances(:, 2) - annotationSelection(i, 3);
                predictionDistances(:, 3) = predictionDistances(:, 3) - annotationSelection(i, 4);
                predictionDistances(:, 4) = predictionDistances(:, 4) - annotationSelection(i, 1);
                predictionDistances = abs(predictionDistances);
                distanceFlags = ...
                    predictionDistances(:, 1) <= searchRadius(1) & ...
                    predictionDistances(:, 2) <= searchRadius(2) & ...
                    predictionDistances(:, 3) <= searchRadius(3) & ...
                    predictionDistances(:, 4) <= searchRadius(4);
                annotationSelection(i, 5) = sum(distanceFlags) > 0;
            end;
            
            annotatedPositives = annotationSelection(:, [2 3 4 1 5]);
            presumedFalseNegatives = annotationSelection(annotationSelection(:, 5) == 0, [2 3 4 1 5]);
            
            save([outputFolder '\Tables\' outputLabel '.Timepoint_' num2str(t) '.annotatedPositives.mat'], 'annotatedPositives');
            save([outputFolder '\Tables\' outputLabel '.Timepoint_' num2str(t) '.presumedFalseNegatives.mat'], 'presumedFalseNegatives');
            
            disp(['extracting image windows around annotations for ' displayLabel ' timepoint ' num2str(t)]);
            
            for y = 1:2
                switch y
                    case 1
                        outputSubfolder = [outputFolder '\' outputLabel '.Timepoint_' num2str(t) '.annotatedPositives'];
                        dataTable = annotatedPositives;
                    case 2
                        outputSubfolder = [outputFolder '\' outputLabel '.Timepoint_' num2str(t) '.presumedFalseNegatives'];
                        dataTable = presumedFalseNegatives;
                end;
                
                if exist(outputSubfolder, 'dir') ~= 7
                    mkdir(outputSubfolder);
                end;
                if saveVolumes && exist([outputSubfolder '\Volumes'], 'dir') ~= 7
                    mkdir([outputSubfolder '\Volumes']);
                end;
                if saveProjections && exist([outputSubfolder '\ProjectionsXY'], 'dir') ~= 7
                    mkdir([outputSubfolder '\ProjectionsXY']);
                end;
                if saveProjections && exist([outputSubfolder '\ProjectionsXZ'], 'dir') ~= 7
                    mkdir([outputSubfolder '\ProjectionsXZ']);
                end;
                
                if ~isempty(dataTable)
                    for i = 1:size(dataTable, 1)
                        minX = max(1, dataTable(i, 1) - croppingRadius(1));
                        maxX = min(size(stack, 1), dataTable(i, 1) + croppingRadius(1));
                        minY = max(1, dataTable(i, 2) - croppingRadius(2));
                        maxY = min(size(stack, 2), dataTable(i, 2) + croppingRadius(2));
                        minZ = max(1, dataTable(i, 3) - croppingRadius(3));
                        maxZ = min(size(stack, 3), dataTable(i, 3) + croppingRadius(3));
                        
                        subVolume = zeros(2*croppingRadius(1)+1, 2*croppingRadius(2)+1, 2*croppingRadius(3)+1, 2*croppingRadius(4)+1, 'uint16');
                        
                        if (dataTable(i, 1) - croppingRadius(1)) < 1
                            xOffset = 1 - (dataTable(i, 1) - croppingRadius(1));
                        else
                            xOffset = 0;
                        end;
                        if (dataTable(i, 2) - croppingRadius(2)) < 1
                            yOffset = 1 - (dataTable(i, 2) - croppingRadius(2));
                        else
                            yOffset = 0;
                        end;
                        if (dataTable(i, 3) - croppingRadius(3)) < 1
                            zOffset = 1 - (dataTable(i, 3) - croppingRadius(3));
                        else
                            zOffset = 0;
                        end;
                        
                        subVolume(...
                            (xOffset + 1):(xOffset + 1 + maxX - minX), ...
                            (yOffset + 1):(yOffset + 1 + maxY - minY), ...
                            (zOffset + 1):(zOffset + 1 + maxZ - minZ), :) = imageBlock(minX:maxX, minY:maxY, minZ:maxZ, :);
                        subProjectionXY = squeeze(max(subVolume, [], 3));
                        subProjectionXZ = squeeze(max(subVolume, [], 2));
                        
                        outputStackName = [outputSubfolder '\Volumes\Volume_' num2str(i, '%.4d') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.klb'];
                        outputProjectionXYName = [outputSubfolder '\ProjectionsXY\ProjectionXY_' num2str(i, '%.4d') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.tif'];
                        outputProjectionXZName = [outputSubfolder '\ProjectionsXZ\ProjectionXZ_' num2str(i, '%.4d') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.tif'];
                        
                        if saveVolumes
                            writeImage(subVolume, outputStackName);
                        end;
                        if saveProjections
                            for tLocal = 1:size(subProjectionXY, 3)
                                if tLocal == 1
                                    writeMode = 'overwrite';
                                else
                                    writeMode = 'append';
                                end;
                                imwrite(subProjectionXY(:, :, tLocal), outputProjectionXYName, 'Compression', 'none', 'WriteMode', writeMode);
                                imwrite(subProjectionXZ(:, :, tLocal), outputProjectionXZName, 'Compression', 'none', 'WriteMode', writeMode);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        
        if processPredictions
            disp(['processing predictions for ' displayLabel ' timepoint ' num2str(t)]);
            
            annotatedPredictionsSelection = annotatedPredictions(annotatedPredictions(:, 4) == t, :);
            annotatedPredictionsSelection = flipud(sortrows(annotatedPredictionsSelection, 5));
            
            for i = 1:size(annotatedPredictionsSelection, 1)
                annotationDistances = annotations;
                annotationDistances(:, 1) = annotationDistances(:, 1) - annotatedPredictionsSelection(i, 4);
                annotationDistances(:, 2) = annotationDistances(:, 2) - annotatedPredictionsSelection(i, 1);
                annotationDistances(:, 3) = annotationDistances(:, 3) - annotatedPredictionsSelection(i, 2);
                annotationDistances(:, 4) = annotationDistances(:, 4) - annotatedPredictionsSelection(i, 3);
                annotationDistances = abs(annotationDistances);
                distanceFlags = ...
                    annotationDistances(:, 1) <= searchRadius(4) & ...
                    annotationDistances(:, 2) <= searchRadius(1) & ...
                    annotationDistances(:, 3) <= searchRadius(2) & ...
                    annotationDistances(:, 4) <= searchRadius(3);
                annotatedPredictionsSelection(i, 6) = sum(distanceFlags) > 0;
            end;
            annotatedPredictionsSelection(:, 1:4) = round(annotatedPredictionsSelection(:, 1:4));
            
            confirmedTruePositives = annotatedPredictionsSelection(annotatedPredictionsSelection(:, 6) == 1, :);
            presumedFalsePositives = annotatedPredictionsSelection(annotatedPredictionsSelection(:, 6) == 0, :);
            
            save([outputFolder '\Tables\' outputLabel '.Timepoint_' num2str(t) '.confirmedTruePositives.mat'], 'confirmedTruePositives');
            save([outputFolder '\Tables\' outputLabel '.Timepoint_' num2str(t) '.presumedFalsePositives.mat'], 'presumedFalsePositives');
            
            disp(['extracting image windows around predictions for ' displayLabel ' timepoint ' num2str(t)]);
            
            for y = 1:2
                switch y
                    case 1
                        outputSubfolder = [outputFolder '\' outputLabel '.Timepoint_' num2str(t) '.confirmedTruePositives'];
                        dataTable = confirmedTruePositives;
                    case 2
                        outputSubfolder = [outputFolder '\' outputLabel '.Timepoint_' num2str(t) '.presumedFalsePositives'];
                        dataTable = presumedFalsePositives;
                end;
                
                if exist(outputSubfolder, 'dir') ~= 7
                    mkdir(outputSubfolder);
                end;
                if saveVolumes && exist([outputSubfolder '\Volumes'], 'dir') ~= 7
                    mkdir([outputSubfolder '\Volumes']);
                end;
                if saveProjections && exist([outputSubfolder '\ProjectionsXY'], 'dir') ~= 7
                    mkdir([outputSubfolder '\ProjectionsXY']);
                end;
                if saveProjections && exist([outputSubfolder '\ProjectionsXZ'], 'dir') ~= 7
                    mkdir([outputSubfolder '\ProjectionsXZ']);
                end;
                
                if ~isempty(dataTable)
                    for i = 1:size(dataTable, 1)
                        minX = max(1, dataTable(i, 1) - croppingRadius(1));
                        maxX = min(size(stack, 1), dataTable(i, 1) + croppingRadius(1));
                        minY = max(1, dataTable(i, 2) - croppingRadius(2));
                        maxY = min(size(stack, 2), dataTable(i, 2) + croppingRadius(2));
                        minZ = max(1, dataTable(i, 3) - croppingRadius(3));
                        maxZ = min(size(stack, 3), dataTable(i, 3) + croppingRadius(3));
                        
                        subVolume = zeros(2*croppingRadius(1)+1, 2*croppingRadius(2)+1, 2*croppingRadius(3)+1, 2*croppingRadius(4)+1, 'uint16');
                        
                        if (dataTable(i, 1) - croppingRadius(1)) < 1
                            xOffset = 1 - (dataTable(i, 1) - croppingRadius(1));
                        else
                            xOffset = 0;
                        end;
                        if (dataTable(i, 2) - croppingRadius(2)) < 1
                            yOffset = 1 - (dataTable(i, 2) - croppingRadius(2));
                        else
                            yOffset = 0;
                        end;
                        if (dataTable(i, 3) - croppingRadius(3)) < 1
                            zOffset = 1 - (dataTable(i, 3) - croppingRadius(3));
                        else
                            zOffset = 0;
                        end;
                        
                        subVolume(...
                            (xOffset + 1):(xOffset + 1 + maxX - minX), ...
                            (yOffset + 1):(yOffset + 1 + maxY - minY), ...
                            (zOffset + 1):(zOffset + 1 + maxZ - minZ), :) = imageBlock(minX:maxX, minY:maxY, minZ:maxZ, :);
                        subProjectionXY = squeeze(max(subVolume, [], 3));
                        subProjectionXZ = squeeze(max(subVolume, [], 2));
                        
                        outputStackName = [outputSubfolder '\Volumes\Volume_' num2str(i, '%.4d') '.Score_' num2str(dataTable(i, 5), '%.8f') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.klb'];
                        outputProjectionXYName = [outputSubfolder '\ProjectionsXY\ProjectionXY_' num2str(i, '%.4d') '.Score_' num2str(dataTable(i, 5), '%.8f') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.tif'];
                        outputProjectionXZName = [outputSubfolder '\ProjectionsXZ\ProjectionXZ_' num2str(i, '%.4d') '.Score_' num2str(dataTable(i, 5), '%.8f') ...
                            '.X' num2str(dataTable(i, 1)) '_Y' num2str(dataTable(i, 2)) '_Z' num2str(dataTable(i, 3)) '.tif'];
                        
                        if saveVolumes
                            writeImage(subVolume, outputStackName);
                        end;
                        if saveProjections
                            for tLocal = 1:size(subProjectionXY, 3)
                                if tLocal == 1
                                    writeMode = 'overwrite';
                                else
                                    writeMode = 'append';
                                end;
                                imwrite(subProjectionXY(:, :, tLocal), outputProjectionXYName, 'Compression', 'none', 'WriteMode', writeMode);
                                imwrite(subProjectionXZ(:, :, tLocal), outputProjectionXZName, 'Compression', 'none', 'WriteMode', writeMode);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        
        disp(' ');
    end;
end;