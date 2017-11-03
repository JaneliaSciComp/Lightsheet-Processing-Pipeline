imageSourceFolder  = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Time-Lapse Image Stacks';
divisionSourceFile = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Annotations\divisionAnnotations.mat';

divisionStacksTimeSeriesOutputFolder      = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Stacks Time Series';
divisionStacksSnapshotsOutputFolder       = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Stacks Snapshots';
divisionProjectionsTimeSeriesOutputFolder = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Projections Time Series';
divisionProjectionsSnapshotsOutputFolder  = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Projections Snapshots';
divisionOverviewOutputFolder              = 'E:' filesep 'Mouse Development' filesep 'Divisions' filesep 'Division Annotations Visualization';

specimen   = 0;
timepoints = 0:531;
cameras    = [0 1];
channels   = 0;
footer     = '.fusedStack.corrected.shifted.klb';

croppingRadii = [40 40 10 3]; % [x, y, z, t]

labelKeys = {...
    1, 'A';...
    2, 'B';...
    3, 'C';...
    4, 'D';...
    5, 'E';...
    0, 'W';...
    103, 'CT';...
    100, 'WT'};

%% data extraction

stackNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM'];
stackNameFooter = '';
for c = cameras
    stackNameFooter = [stackNameFooter 'CM' num2str(c, '%.2d') '_'];
end;
for c = channels
    if c == channels(end)
        stackNameFooter = [stackNameFooter 'CHN' num2str(c, '%.2d')];
    else
        stackNameFooter = [stackNameFooter 'CHN' num2str(c, '%.2d') '_'];
    end;
end;

if exist(divisionStacksTimeSeriesOutputFolder, 'dir') ~= 7
    mkdir(divisionStacksTimeSeriesOutputFolder);
end;
if exist(divisionStacksSnapshotsOutputFolder, 'dir') ~= 7
    mkdir(divisionStacksSnapshotsOutputFolder);
end;
if exist(divisionProjectionsTimeSeriesOutputFolder, 'dir') ~= 7
    mkdir(divisionProjectionsTimeSeriesOutputFolder);
end;
if exist(divisionProjectionsSnapshotsOutputFolder, 'dir') ~= 7
    mkdir(divisionProjectionsSnapshotsOutputFolder);
end;
if exist(divisionOverviewOutputFolder, 'dir') ~= 7
    mkdir(divisionOverviewOutputFolder);
end;

load(divisionSourceFile);

divisionAnnotations = sortrows(divisionAnnotations(:, 1:5), 1);
divisionAnnotations(:, [3 4 5]) = round(divisionAnnotations(:, [3 4 5]));

minT = divisionAnnotations(1, 1);
maxT = divisionAnnotations(end, 1);
uniqueTimepoints = unique(divisionAnnotations(:, 1));

annotationLabels = cell(size(divisionAnnotations, 1), 1);
labelKeysNumerical = cell2mat(labelKeys(:, 1));
for i = 1:size(divisionAnnotations, 1)
    annotationLabels{i, 1} = labelKeys{find(labelKeysNumerical == divisionAnnotations(i, 2), 1), 2};
end;

croppedVolumes = zeros(croppingRadii(1)*2+1, croppingRadii(2)*2+1, croppingRadii(3)*2+1, croppingRadii(4)*2+1, size(divisionAnnotations, 1), 'uint16');
volumeSlots = ones(size(divisionAnnotations, 1), 1);

croppedXYProjections = zeros(croppingRadii(1)*2+1, croppingRadii(2)*2+1, croppingRadii(4)*2+1, size(divisionAnnotations, 1), 'uint16');
croppedXZProjections = zeros(croppingRadii(1)*2+1, croppingRadii(3)*2+1, croppingRadii(4)*2+1, size(divisionAnnotations, 1), 'uint16');
croppedYZProjections = zeros(croppingRadii(2)*2+1, croppingRadii(3)*2+1, croppingRadii(4)*2+1, size(divisionAnnotations, 1), 'uint16');

overviewXYProjections = zeros(croppingRadii(1)*2+1, croppingRadii(2)*2+1, size(divisionAnnotations, 1), 'uint16');
overviewXZProjections = zeros(croppingRadii(1)*2+1, croppingRadii(3)*2+1, size(divisionAnnotations, 1), 'uint16');
overviewYZProjections = zeros(croppingRadii(2)*2+1, croppingRadii(3)*2+1, size(divisionAnnotations, 1), 'uint16');

tic;

disp(' ');
for t = timepoints
    proximityFlag = sum(abs(uniqueTimepoints - t) <= croppingRadii(4)) > 0;
    if proximityFlag
        sourceStackName = [stackNameHeader num2str(t, '%.6d') '_' stackNameFooter footer];
        sourceStackFullName = [imageSourceFolder filesep '' sourceStackName];
        
        disp(['Loading and processing source volume at time point ' num2str(t)]);
        currentStack = readImage(sourceStackFullName);
        [xSize, ySize, zSize] = size(currentStack);
        
        currentIndices = find(abs(divisionAnnotations(:, 1) - t) <= croppingRadii(4));
        for i = currentIndices'
            currentXCenter = divisionAnnotations(i, 3);
            currentYCenter = divisionAnnotations(i, 4);
            currentZCenter = divisionAnnotations(i, 5);
            currentSubvolume = currentStack(...
                max(1, (currentXCenter-croppingRadii(1))):min(xSize, (currentXCenter+croppingRadii(1))),...
                max(1, (currentYCenter-croppingRadii(2))):min(ySize, (currentYCenter+croppingRadii(2))),...
                max(1, (currentZCenter-croppingRadii(3))):min(zSize, (currentZCenter+croppingRadii(3))));
            [xSubvolumeSize, ySubvolumeSize, zSubvolumeSize] = size(currentSubvolume);
            if (currentXCenter - croppingRadii(1)) < 1
                lowXPadding = croppingRadii(1) - currentXCenter + 1;
                currentSubvolume = cat(1, zeros(lowXPadding, ySubvolumeSize, zSubvolumeSize, 'uint16'), currentSubvolume);
                xSubvolumeSize = xSubvolumeSize + lowXPadding;
            end;
            if (currentXCenter + croppingRadii(1)) > xSize
                highXPadding = currentXCenter + croppingRadii(1) - xSize;
                currentSubvolume = cat(1, currentSubvolume, zeros(highXPadding, ySubvolumeSize, zSubvolumeSize, 'uint16'));
                xSubvolumeSize = xSubvolumeSize + highXPadding;
            end;
            if (currentYCenter - croppingRadii(2)) < 1
                lowYPadding = croppingRadii(2) - currentYCenter + 1;
                currentSubvolume = cat(2, zeros(xSubvolumeSize, lowYPadding, zSubvolumeSize, 'uint16'), currentSubvolume);
                ySubvolumeSize = ySubvolumeSize + lowYPadding;
            end;
            if (currentYCenter + croppingRadii(2)) > ySize
                highYPadding = currentYCenter + croppingRadii(2) - ySize;
                currentSubvolume = cat(2, currentSubvolume, zeros(xSubvolumeSize, highYPadding, zSubvolumeSize, 'uint16'));
                ySubvolumeSize = ySubvolumeSize + highYPadding;
            end;
            if (currentZCenter - croppingRadii(3)) < 1
                lowZPadding = croppingRadii(3) - currentZCenter + 1;
                currentSubvolume = cat(3, zeros(xSubvolumeSize, ySubvolumeSize, lowZPadding, 'uint16'), currentSubvolume);
                zSubvolumeSize = zSubvolumeSize + lowZPadding;
            end;
            if (currentZCenter + croppingRadii(3)) > zSize
                highZPadding = currentZCenter + croppingRadii(3) - zSize;
                currentSubvolume = cat(3, currentSubvolume, zeros(xSubvolumeSize, ySubvolumeSize, highZPadding, 'uint16'));
                zSubvolumeSize = zSubvolumeSize + highZPadding;
            end;
            croppedVolumes(:, :, :, volumeSlots(i), i) = currentSubvolume;
            
            croppedXYProjections(:, :, volumeSlots(i), i) = max(currentSubvolume, [], 3);
            croppedXZProjections(:, :, volumeSlots(i), i) = squeeze(max(currentSubvolume, [], 2));
            croppedYZProjections(:, :, volumeSlots(i), i) = squeeze(max(currentSubvolume, [], 1));
            
            if divisionAnnotations(i, 1) == t
                currentDivisionStackSnapshotFullName = [divisionStacksSnapshotsOutputFolder ...
                    filesep 'StackSnapshot.Class_' annotationLabels{i, 1} '.TM_' num2str(t, '%.4d') ...
                    '.Centroid_' num2str(currentXCenter, '%.4d') '_' num2str(currentYCenter, '%.4d') '_' num2str(currentZCenter, '%.4d') '.klb'];
                writeImage(currentSubvolume, currentDivisionStackSnapshotFullName);
                
                currentDivisionProjectionXYSnapshotFullName = [divisionProjectionsSnapshotsOutputFolder ...
                    filesep 'ProjectionXYSnapshot.Class_' annotationLabels{i, 1} '.TM_' num2str(t, '%.4d') ...
                    '.Centroid_' num2str(currentXCenter, '%.4d') '_' num2str(currentYCenter, '%.4d') '_' num2str(currentZCenter, '%.4d') '.klb'];
                writeImage(croppedXYProjections(:, :, volumeSlots(i), i), currentDivisionProjectionXYSnapshotFullName);
                
                currentDivisionProjectionXZSnapshotFullName = [divisionProjectionsSnapshotsOutputFolder ...
                    filesep 'ProjectionXZSnapshot.Class_' annotationLabels{i, 1} '.TM_' num2str(t, '%.4d') ...
                    '.Centroid_' num2str(currentXCenter, '%.4d') '_' num2str(currentYCenter, '%.4d') '_' num2str(currentZCenter, '%.4d') '.klb'];
                writeImage(croppedXZProjections(:, :, volumeSlots(i), i), currentDivisionProjectionXZSnapshotFullName);
                
                currentDivisionProjectionYZSnapshotFullName = [divisionProjectionsSnapshotsOutputFolder ...
                    filesep 'ProjectionYZSnapshot.Class_' annotationLabels{i, 1} '.TM_' num2str(t, '%.4d') ...
                    '.Centroid_' num2str(currentXCenter, '%.4d') '_' num2str(currentYCenter, '%.4d') '_' num2str(currentZCenter, '%.4d') '.klb'];
                writeImage(croppedYZProjections(:, :, volumeSlots(i), i), currentDivisionProjectionYZSnapshotFullName);
                
                overviewXYProjections(:, :, i) = croppedXYProjections(:, :, volumeSlots(i), i);
                overviewXZProjections(:, :, i) = croppedXZProjections(:, :, volumeSlots(i), i);
                overviewYZProjections(:, :, i) = croppedYZProjections(:, :, volumeSlots(i), i);
            end;
            
            volumeSlots(i) = volumeSlots(i) + 1;
        end;
    end;
end;

disp(' ');
disp('Writing overview projection stacks to disk');

overviewDivisionProjectionXYSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXYOverviewSnapshots.klb'];
writeImage(overviewXYProjections, overviewDivisionProjectionXYSnapshotsFullName);

overviewDivisionProjectionXZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXZOverviewSnapshots.klb'];
writeImage(overviewXZProjections, overviewDivisionProjectionXZSnapshotsFullName);

overviewDivisionProjectionYZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionYZOverviewSnapshots.klb'];
writeImage(overviewYZProjections, overviewDivisionProjectionYZSnapshotsFullName);

divisionAnnotationsByClass = cell(size(labelKeys, 1), 1);
indicesByClass = cell(size(labelKeys, 1), 1);
for i = 1:size(labelKeys, 1)
    divisionAnnotationsByClass{i, 1} = divisionAnnotations(divisionAnnotations(:, 2) == labelKeys{i, 1}, :);
    indicesByClass{i, 1} = find(divisionAnnotations(:, 2) == labelKeys{i, 1});
    
    currentClassOverviewDivisionProjectionXYSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXYOverviewSnapshots.Class_' labelKeys{i, 2} '.klb'];
    writeImage(overviewXYProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionXYSnapshotsFullName);
    
    currentClassOverviewDivisionProjectionXZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXZOverviewSnapshots.Class_' labelKeys{i, 2} '.klb'];
    writeImage(overviewXZProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionXZSnapshotsFullName);
    
    currentClassOverviewDivisionProjectionYZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionYZOverviewSnapshots.Class_' labelKeys{i, 2} '.klb'];
    writeImage(overviewYZProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionYZSnapshotsFullName);
end;

for i = 1:size(divisionAnnotations, 1)
    overviewXYProjections(:, :, i) = overviewXYProjections(:, :, i) - min(min(overviewXYProjections(:, :, i)));
    overviewXYProjections(:, :, i) = uint16(double(overviewXYProjections(:, :, i)) .* (65535 / double(max(max(overviewXYProjections(:, :, i))))));
    
    overviewXZProjections(:, :, i) = overviewXZProjections(:, :, i) - min(min(overviewXZProjections(:, :, i)));
    overviewXZProjections(:, :, i) = uint16(double(overviewXZProjections(:, :, i)) .* (65535 / double(max(max(overviewXZProjections(:, :, i))))));
    
    overviewYZProjections(:, :, i) = overviewYZProjections(:, :, i) - min(min(overviewYZProjections(:, :, i)));
    overviewYZProjections(:, :, i) = uint16(double(overviewYZProjections(:, :, i)) .* (65535 / double(max(max(overviewYZProjections(:, :, i))))));
end;

overviewDivisionProjectionXYSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXYOverviewSnapshots.normalized.klb'];
writeImage(overviewXYProjections, overviewDivisionProjectionXYSnapshotsFullName);

overviewDivisionProjectionXZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXZOverviewSnapshots.normalized.klb'];
writeImage(overviewXZProjections, overviewDivisionProjectionXZSnapshotsFullName);

overviewDivisionProjectionYZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionYZOverviewSnapshots.normalized.klb'];
writeImage(overviewYZProjections, overviewDivisionProjectionYZSnapshotsFullName);

divisionAnnotationsByClass = cell(size(labelKeys, 1), 1);
indicesByClass = cell(size(labelKeys, 1), 1);
for i = 1:size(labelKeys, 1)
    divisionAnnotationsByClass{i, 1} = divisionAnnotations(divisionAnnotations(:, 2) == labelKeys{i, 1}, :);
    indicesByClass{i, 1} = find(divisionAnnotations(:, 2) == labelKeys{i, 1});
    
    currentClassOverviewDivisionProjectionXYSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXYOverviewSnapshots.Class_' labelKeys{i, 2} '.normalized.klb'];
    writeImage(overviewXYProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionXYSnapshotsFullName);
    
    currentClassOverviewDivisionProjectionXZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionXZOverviewSnapshots.Class_' labelKeys{i, 2} '.normalized.klb'];
    writeImage(overviewXZProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionXZSnapshotsFullName);
    
    currentClassOverviewDivisionProjectionYZSnapshotsFullName = [divisionOverviewOutputFolder filesep 'ProjectionYZOverviewSnapshots.Class_' labelKeys{i, 2} '.normalized.klb'];
    writeImage(overviewYZProjections(:, :, indicesByClass{i, 1}'), currentClassOverviewDivisionProjectionYZSnapshotsFullName);
end;

disp(' ');
disp('Writing time series image data to disk');

for i = 1:size(divisionAnnotations, 1)
    currentDivisionStackTimeSeriesFullName = [divisionStacksTimeSeriesOutputFolder ...
        filesep 'StackTimeSeries.Class_' annotationLabels{i, 1} '.TM_' num2str(divisionAnnotations(i, 1), '%.4d') ...
        '.Centroid_' num2str(divisionAnnotations(i, 3), '%.4d') '_' num2str(divisionAnnotations(i, 4), '%.4d') '_' num2str(divisionAnnotations(i, 5), '%.4d') '.klb'];
    writeImage(croppedVolumes(:, :, :, :, i), currentDivisionStackTimeSeriesFullName);
    
    currentDivisionProjectionXYTimeSeriesFullName = [divisionProjectionsTimeSeriesOutputFolder ...
        filesep 'ProjectionXYTimeSeries.Class_' annotationLabels{i, 1} '.TM_' num2str(divisionAnnotations(i, 1), '%.4d') ...
        '.Centroid_' num2str(divisionAnnotations(i, 3), '%.4d') '_' num2str(divisionAnnotations(i, 4), '%.4d') '_' num2str(divisionAnnotations(i, 5), '%.4d') '.klb'];
    writeImage(croppedXYProjections(:, :, :, i), currentDivisionProjectionXYTimeSeriesFullName);
    
    currentDivisionProjectionXZTimeSeriesFullName = [divisionProjectionsTimeSeriesOutputFolder ...
        filesep 'ProjectionXZTimeSeries.Class_' annotationLabels{i, 1} '.TM_' num2str(divisionAnnotations(i, 1), '%.4d') ...
        '.Centroid_' num2str(divisionAnnotations(i, 3), '%.4d') '_' num2str(divisionAnnotations(i, 4), '%.4d') '_' num2str(divisionAnnotations(i, 5), '%.4d') '.klb'];
    writeImage(croppedXZProjections(:, :, :, i), currentDivisionProjectionXZTimeSeriesFullName);
    
    currentDivisionProjectionYZTimeSeriesFullName = [divisionProjectionsTimeSeriesOutputFolder ...
        filesep 'ProjectionYZTimeSeries.Class_' annotationLabels{i, 1} '.TM_' num2str(divisionAnnotations(i, 1), '%.4d') ...
        '.Centroid_' num2str(divisionAnnotations(i, 3), '%.4d') '_' num2str(divisionAnnotations(i, 4), '%.4d') '_' num2str(divisionAnnotations(i, 5), '%.4d') '.klb'];
    writeImage(croppedYZProjections(:, :, :, i), currentDivisionProjectionYZTimeSeriesFullName);
end;

disp(' ');

volumeCounts = volumeSlots - 1;
save('processingLog.mat', 'divisionAnnotations', 'annotationLabels', 'volumeCounts', 'divisionAnnotationsByClass', 'indicesByClass', 'labelKeys');

elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' s']);
disp(' ');