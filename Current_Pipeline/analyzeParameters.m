function analyzeParameters(...
    timepoints, fullInterval,...
    outputString, outputID, dataType,...
    specimen, cameras, channels, readFactors, smoothing,...
    offsetRange, angleRange, intRange, averaging, staticFlag, configRoot)

% ----------------------------------------------------------------------
% | Time-lapse image data fusion parameter analysis for timeFuse.m     |
% |                                                                    |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2015 |
% | Email: kellerp@janelia.hhmi.org                                    |
% ----------------------------------------------------------------------

% lookUpTable format (4-view fusion): (will be reduced to rows 2-14)
% row 1: time point
% row 2: camera 0 bestOffset
% row 3: camera 0 bestRotation
% row 4: camera 0 correctionFactor
% row 5: camera 0 correctionFlag
% row 6: camera 1 bestOffset
% row 7: camera 1 bestRotation
% row 8: camera 1 correctionFactor
% row 9: camera 1 correctionFlag
% row 10: final bestXOffset
% row 11: final bestYOffset
% row 12: final bestRotation
% row 13: final correctionFactor
% row 14: final correctionFlag

% lookUpTable format (2-view channel fusion): (will be reduced to rows 2-5)
% row 1: time point
% row 2: camera X bestOffset
% row 3: camera X bestRotation
% row 4: camera X correctionFactor
% row 5: camera X correctionFlag

% lookUpTable format (2-view camera fusion): (will be reduced to rows 2-6)
% row 1: time point
% row 2: final bestXOffset
% row 3: final bestYOffset
% row 4: final bestRotation
% row 5: final correctionFactor
% row 6: final correctionFlag

if nargin<16
    configRoot = '';
end

if length(cameras) == 2 && length(channels) == 2
    processingMode = 0; % 4-view fusion
    transformationEntries = 8;
    intensityEntries = 7;
    lookUpTableEntries = 14;
    globalHeader = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')];
elseif length(cameras) == 1 && length(channels) == 2
    processingMode = 1; % 2-view channel fusion
    transformationEntries = 3;
    intensityEntries = 3;
    lookUpTableEntries = 5;
    globalHeader = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')];
elseif length(cameras) == 2 && length(channels) == 1
    processingMode = 2; % 2-view camera fusion
    transformationEntries = 4;
    intensityEntries = 3;
    lookUpTableEntries = 6;
    globalHeader = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')];
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion)');
end;

transformations = NaN(length(timepoints), transformationEntries);
transformations(:, 1) = timepoints';
if readFactors
    intensityCorrections = NaN(length(timepoints), intensityEntries);
    intensityCorrections(:, 1) = timepoints';
end;
lookUpTable = NaN(length(fullInterval), lookUpTableEntries);
lookUpTable(:, 1) = fullInterval';

resultsDirectory = [configRoot filesep 'SPM' num2str(specimen, '%.2d') globalHeader '_analyzeParameters'];
mkdir(resultsDirectory);

for currentTP = timepoints
    outputFolder = [outputString '.TM' num2str(currentTP, '%.6d') '_multiFused' outputID];
    currentIndex = find(timepoints == currentTP, 1);
    
    switch processingMode
        case 0            
            for camera = cameras
                transformationIndex = 2 + 2 * (find(cameras == camera, 1) - 1);
                transformationName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') '.transformation.mat'];
                if exist(transformationName, 'file') == 2
                    load(transformationName);
                    transformations(currentIndex, transformationIndex:(transformationIndex + 1)) = transformation(2:3);
                end;
            end;
            
            transformationName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') '.transformation.mat'];
            if exist(transformationName, 'file') == 2
                load(transformationName);
                transformations(currentIndex, 6:8) = transformation;
            end;
            
            if readFactors
                for camera = cameras
                    intensityCorrectionIndex = 2 + 2 * (find(cameras == camera, 1) - 1);
                    intensityCorrectionName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') '_CM' num2str(camera, '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') '.intensityCorrection.mat'];
                    if exist(intensityCorrectionName, 'file') == 2
                        load(intensityCorrectionName);
                        if dataType == 0
                            intensityCorrections(currentIndex, intensityCorrectionIndex:(intensityCorrectionIndex + 1)) = intensityCorrection(4:5);
                        else
                            intensityCorrections(currentIndex, intensityCorrectionIndex:(intensityCorrectionIndex + 1)) = intensityCorrection(5:6);
                        end;
                    end;
                end;
                
                intensityCorrectionName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') '_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') '.intensityCorrection.mat'];
                if exist(intensityCorrectionName, 'file') == 2
                    load(intensityCorrectionName);
                    intensityCorrections(currentIndex, 6:7) = intensityCorrection(4:5);
                end;
            end;
        case 1
            transformationName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') globalHeader '.transformation.mat'];
            if exist(transformationName, 'file') == 2
                load(transformationName);
                transformations(currentIndex, 2:3) = transformation(2:3);
            end;
            
            intensityCorrectionName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') globalHeader '.intensityCorrection.mat'];
            if exist(intensityCorrectionName, 'file') == 2
                load(intensityCorrectionName);
                if dataType == 0
                    intensityCorrections(currentIndex, 2:3) = intensityCorrection(4:5);
                else
                    intensityCorrections(currentIndex, 2:3) = intensityCorrection(5:6);
                end;
            end;
        case 2
            transformationName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') globalHeader '.transformation.mat'];
            if exist(transformationName, 'file') == 2
                load(transformationName);
                transformations(currentIndex, 2:4) = transformation;
            end;
            
            intensityCorrectionName = [outputFolder filesep 'SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') globalHeader '.intensityCorrection.mat'];
            if exist(intensityCorrectionName, 'file') == 2
                load(intensityCorrectionName);
                if dataType == 0
                    intensityCorrections(currentIndex, 2:3) = intensityCorrection(4:5);
                else
                    intensityCorrections(currentIndex, 2:3) = intensityCorrection(5:6);
                end;
            end;
    end;
end;

save([resultsDirectory filesep 'transformations.mat'], 'transformations');

if readFactors
    intensityCorrections(intensityCorrections(:, 3) == 1, 2) = 1 ./ intensityCorrections(intensityCorrections(:, 3) == 1, 2);
    intensityCorrections(intensityCorrections(:, 3) == 1, 3) = 2;
    
    if processingMode == 0
        intensityCorrections(intensityCorrections(:, 5) == 1, 4) = 1 ./ intensityCorrections(intensityCorrections(:, 5) == 1, 4);
        intensityCorrections(intensityCorrections(:, 5) == 1, 5) = 2;
        
        intensityCorrections(intensityCorrections(:, 7) == 1, 6) = 1 ./ intensityCorrections(intensityCorrections(:, 7) == 1, 6);
        intensityCorrections(intensityCorrections(:, 7) == 1, 7) = 2;
    end;
    
    save([resultsDirectory filesep 'intensityCorrections.mat'], 'intensityCorrections');
end;

% visualize parameters, create median-filtered and mean-filtered lookUpTable

fh=figure('visible','off');
h=axes();
switch processingMode
    case 0
        plot(transformations(:, 1), transformations(:, 2), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        hold on;
        plot(transformations(:, 1), transformations(:, 4), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        plot(transformations(:, 1), transformations(:, 6), 'Color', [0.5 0 1], 'LineWidth', 1.5);
        plot(transformations(:, 1), transformations(:, 7), 'Color', [0 0.5 1], 'LineWidth', 1.5);
    case 1
        if cameras == 0
            plot(transformations(:, 1), transformations(:, 2), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        else
            plot(transformations(:, 1), transformations(:, 2), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        end;
    case 2
        plot(transformations(:, 1), transformations(:, 2), 'Color', [0.5 0 1], 'LineWidth', 1.5);
        plot(transformations(:, 1), transformations(:, 3), 'Color', [0 0.5 1], 'LineWidth', 1.5);
end;
set(h, 'FontSize', 14, 'LineWidth', 1.5);
xlabel(h, 'time point'); ylabel(h, 'offset (pixels)');
currentFile = [resultsDirectory filesep 'offsets.png']; saveas(h, currentFile);
currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
close(fh);

fh=figure('visible','off');
h=axes();
switch processingMode
    case 0
        plot(transformations(:, 1), transformations(:, 3), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        hold on;
        plot(transformations(:, 1), transformations(:, 5), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        plot(transformations(:, 1), transformations(:, 8), 'Color', [0 0 1], 'LineWidth', 1.5);
    case 1
        if cameras == 0
            plot(transformations(:, 1), transformations(:, 3), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        else
            plot(transformations(:, 1), transformations(:, 3), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        end;
    case 2
        plot(transformations(:, 1), transformations(:, 4), 'Color', [0 0 1], 'LineWidth', 1.5);
end;
set(h, 'FontSize', 14, 'LineWidth', 1.5);
xlabel(h, 'time point'); ylabel(h, 'angle (degrees)');
currentFile = [resultsDirectory filesep 'angles.png']; saveas(h, currentFile);
currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
close(fh);

fh=figure('visible','off');
h=axes();
switch processingMode
    case 0
        plot(intensityCorrections(:, 1), intensityCorrections(:, 2), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        hold on;
        plot(intensityCorrections(:, 1), intensityCorrections(:, 4), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        plot(intensityCorrections(:, 1), intensityCorrections(:, 6), 'Color', [0 0 1], 'LineWidth', 1.5);
    case 1
        if cameras == 0
            plot(intensityCorrections(:, 1), intensityCorrections(:, 2), 'Color', [1 0.5 0], 'LineWidth', 1.5);
        else
            plot(intensityCorrections(:, 1), intensityCorrections(:, 2), 'Color', [1 0 0.5], 'LineWidth', 1.5);
        end;
    case 2
        plot(intensityCorrections(:, 1), intensityCorrections(:, 2), 'Color', [0 0 1], 'LineWidth', 1.5);
end;
set(h, 'FontSize', 14, 'LineWidth', 1.5);
xlabel(h, 'time point'); ylabel(h, 'factor');
currentFile = [resultsDirectory filesep 'factors.png']; saveas(h, currentFile);
currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
close(fh);

% load transformations; load intensityCorrections;

smoothTransformations = transformations;
smoothIntensityCorrections = intensityCorrections;

smoothingRange1 = 2:transformationEntries;
smoothingRange2 = 2:2:intensityEntries;

switch processingMode
    case 0
        lookUpTableSlots1 = [2:3 6:7 10:12];
        lookUpTableSlots2 = [4 8 13];
    case 1
        lookUpTableSlots1 = 2:3;
        lookUpTableSlots2 = 4;
    case 2
        lookUpTableSlots1 = 2:4;
        lookUpTableSlots2 = 5;
end;

if smoothing(1)
    for i = smoothingRange1
        smoothTransformations(:, i) = smooth(transformations(:, i), smoothing(2), 'rloess');
    end;
    if readFactors
        for i = smoothingRange2
            smoothIntensityCorrections(:, i) = smooth(intensityCorrections(:, i), smoothing(2), 'rloess');
        end;
    end;
end;

for currentTP = fullInterval
    currentIndex = find(timepoints == currentTP, 1);
    if ~isempty(currentIndex)
        lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots1) = smoothTransformations(currentIndex, smoothingRange1);
        if readFactors
            lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots2) = smoothIntensityCorrections(currentIndex, smoothingRange2);
        end;
    else
        distances = timepoints - currentTP;
        positiveDistances = distances .* (distances > 0);
        negativeDistances = distances .* (distances < 0);
        if ~isempty(positiveDistances(positiveDistances > 0)) && ~isempty(negativeDistances(negativeDistances < 0))
            upperIndex = find(positiveDistances == min(positiveDistances(positiveDistances > 0)), 1);
            lowerIndex = find(negativeDistances == max(negativeDistances(negativeDistances < 0)), 1);
            upperDistance = min(positiveDistances(positiveDistances > 0));
            lowerDistance = -max(negativeDistances(negativeDistances < 0));
            upperWeighting = lowerDistance / (lowerDistance + upperDistance);
            lowerWeighting = upperDistance / (lowerDistance + upperDistance);
            lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots1) = upperWeighting .* smoothTransformations(upperIndex, smoothingRange1) + lowerWeighting .* smoothTransformations(lowerIndex, smoothingRange1);
            if readFactors
                lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots2) = upperWeighting .* smoothIntensityCorrections(upperIndex, smoothingRange2) + lowerWeighting .* smoothIntensityCorrections(lowerIndex, smoothingRange2);
            end;
        else
            distances = abs(distances);
            closestIndex = find(distances == min(distances), 1);
            lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots1) = smoothTransformations(closestIndex, smoothingRange1);
            if readFactors
                lookUpTable(find(fullInterval == currentTP, 1), lookUpTableSlots2) = smoothIntensityCorrections(closestIndex, smoothingRange2);
            end;
        end;
    end;
end;

medianTable = NaN(length(fullInterval), lookUpTableEntries);
medianTable(:, 1) = fullInterval';
meanTable = NaN(length(fullInterval), lookUpTableEntries);
meanTable(:, 1) = fullInterval';

switch processingMode
    case 0
        lookUpTableSlots1 = [2 6 10 11];
        lookUpTableSlots2 = [3 7 12];
        lookUpTableSlots3 = [4 8 13];
        lookUpTableSlots4 = [5 9 14];
    case 1
        lookUpTableSlots1 = 2;
        lookUpTableSlots2 = 3;
        lookUpTableSlots3 = 4;
        lookUpTableSlots4 = 5;
    case 2
        lookUpTableSlots1 = 2:3;
        lookUpTableSlots2 = 4;
        lookUpTableSlots3 = 5;
        lookUpTableSlots4 = 6;
end;

for currentIndex = 1:length(fullInterval)
    currentOffsetRange = max(1, currentIndex - offsetRange):min(length(fullInterval), currentIndex + offsetRange);
    medianTable(currentIndex, lookUpTableSlots1) = median(lookUpTable(currentOffsetRange, lookUpTableSlots1), 1);
    meanTable(currentIndex, lookUpTableSlots1)   = mean(lookUpTable(currentOffsetRange, lookUpTableSlots1), 1);
    
    currentAngleRange = max(1, currentIndex - angleRange):min(length(fullInterval), currentIndex + angleRange);
    medianTable(currentIndex, lookUpTableSlots2) = median(lookUpTable(currentAngleRange, lookUpTableSlots2), 1);
    meanTable(currentIndex, lookUpTableSlots2)   = mean(lookUpTable(currentAngleRange, lookUpTableSlots2), 1);
    
    currentIntRange = max(1, currentIndex - intRange):min(length(fullInterval), currentIndex + intRange);
    medianTable(currentIndex, lookUpTableSlots3) = median(lookUpTable(currentIntRange, lookUpTableSlots3), 1);
    meanTable(currentIndex, lookUpTableSlots3)   = mean(lookUpTable(currentIntRange, lookUpTableSlots3), 1);
    medianTable(currentIndex, lookUpTableSlots4) = lookUpTable(currentIndex, lookUpTableSlots4);
    meanTable(currentIndex, lookUpTableSlots4)   = lookUpTable(currentIndex, lookUpTableSlots4);
end;

if staticFlag
    staticMedian = mean(medianTable(:, 2:end), 1);
    staticMean = mean(meanTable(:, 2:end), 1);
    for currentIndex = 1:length(fullInterval)
        medianTable(currentIndex, 2:end) = staticMedian;
        meanTable(currentIndex, 2:end) = staticMean;
    end;
end;

switch processingMode
    case 0
        if readFactors
            indices = [2:4 6:8 10:13];
        else
            indices = [2:3 6:7 10:12];
        end;
        
        for currentIndex = indices
            switch currentIndex
                case 2
                    header = ['CM' num2str(cameras(1)) '_bestOffset'];
                    label  = ['CM' num2str(cameras(1)) ' bestOffset'];
                    unit   = ' (pixels)';
                    oldIndex = 2;
                case 3
                    header = ['CM' num2str(cameras(1)) '_bestRotation'];
                    label  = ['CM' num2str(cameras(1)) 'bestRotation'];
                    unit   = ' (degrees)';
                    oldIndex = 3;
                case 4
                    header = ['CM' num2str(cameras(1)) '_correctionFactor'];
                    label  = ['CM' num2str(cameras(1)) 'correctionFactor'];
                    unit   = '';
                    oldIndex = 2;
                case 6
                    header = ['CM' num2str(cameras(2)) '_bestOffset'];
                    label  = ['CM' num2str(cameras(2)) ' bestOffset'];
                    unit   = ' (pixels)';
                    oldIndex = 4;
                case 7
                    header = ['CM' num2str(cameras(2)) '_bestRotation'];
                    label  = ['CM' num2str(cameras(2)) 'bestRotation'];
                    unit   = ' (degrees)';
                    oldIndex = 5;
                case 8
                    header = ['CM' num2str(cameras(2)) '_correctionFactor'];
                    label  = ['CM' num2str(cameras(2)) 'correctionFactor'];
                    unit   = '';
                    oldIndex = 4;
                case 10
                    header = 'final_bestXOffset';
                    label  = 'final bestXOffset';
                    unit   = ' (pixels)';
                    oldIndex = 6;
                case 11
                    header = 'final_bestYOffset';
                    label  = 'final bestYOffset';
                    unit   = ' (pixels)';
                    oldIndex = 7;
                case 12
                    header = 'final_bestRotation';
                    label  = 'final bestRotation';
                    unit   = ' (degrees)';
                    oldIndex = 8;
                case 13
                    header = 'final_correctionFactor';
                    label  = 'final correctionFactor';
                    unit   = '';
                    oldIndex = 6;
            end;
            
            fh=figure('visible','off');
            h=axes();
            if ~isempty(find([2:3 6:7 10:12] == currentIndex, 1))
                plot(transformations(:, 1), transformations(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            else
                plot(intensityCorrections(:, 1), intensityCorrections(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            end;
            hold on;
            plot(medianTable(:, 1), medianTable(:, currentIndex), 'Color', [1 0 0], 'LineWidth', 1.5);
            plot(meanTable(:, 1), meanTable(:, currentIndex), '--', 'Color', [1 0 0], 'LineWidth', 1.5);
            set(h, 'FontSize', 14, 'LineWidth', 1.5);
            xlabel(h, 'time point'); ylabel(h, [label unit]);
            currentFile = [resultsDirectory filesep '' header '.png']; saveas(h, currentFile);
            currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
            close(fh);
            % close;
        end;
    case 1
        if readFactors
            indices = 2:4;
        else
            indices = 2:3;
        end;
        
        for currentIndex = indices
            switch currentIndex
                case 2
                    header = ['CM' num2str(cameras) '_bestOffset'];
                    label  = ['CM' num2str(cameras) ' bestOffset'];
                    unit   = ' (pixels)';
                    oldIndex = 2;
                case 3
                    header = ['CM' num2str(cameras) '_bestRotation'];
                    label  = ['CM' num2str(cameras) ' bestRotation'];
                    unit   = ' (degrees)';
                    oldIndex = 3;
                case 4
                    header = ['CM' num2str(cameras) '_correctionFactor'];
                    label  = ['CM' num2str(cameras) ' correctionFactor'];
                    unit   = '';
                    oldIndex = 2;
            end;
            
            fh=figure('visible','off');
            h=axes();
            if ~isempty(find(2:3 == currentIndex, 1))
                plot(transformations(:, 1), transformations(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            else
                plot(intensityCorrections(:, 1), intensityCorrections(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            end;
            hold on;
            plot(medianTable(:, 1), medianTable(:, currentIndex), 'Color', [1 0 0], 'LineWidth', 1.5);
            plot(meanTable(:, 1), meanTable(:, currentIndex), '--', 'Color', [1 0 0], 'LineWidth', 1.5);
            set(h, 'FontSize', 14, 'LineWidth', 1.5);
            xlabel(h, 'time point'); ylabel(h, [label unit]);
            currentFile = [resultsDirectory filesep '' header '.png']; saveas(h, currentFile);
            currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
            close(fh);
            % close;
        end;
    case 2
        if readFactors
            indices = 2:5;
        else
            indices = 2:4;
        end;
        
        for currentIndex = indices
            switch currentIndex
                case 2
                    header = 'final_bestXOffset';
                    label  = 'final bestXOffset';
                    unit   = ' (pixels)';
                    oldIndex = 2;
                case 3
                    header = 'final_bestYOffset';
                    label  = 'final bestYOffset';
                    unit   = ' (pixels)';
                    oldIndex = 3;
                case 4
                    header = 'final_bestRotation';
                    label  = 'final bestRotation';
                    unit   = ' (degrees)';
                    oldIndex = 4;
                case 5
                    header = 'final_correctionFactor';
                    label  = 'final correctionFactor';
                    unit   = '';
                    oldIndex = 2;
            end;
            
            fh=figure('visible','off');
            h=axes();
            if ~isempty(find(2:4 == currentIndex, 1))
                plot(transformations(:, 1), transformations(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            else
                plot(intensityCorrections(:, 1), intensityCorrections(:, oldIndex), 'Color', [0 0 1], 'LineWidth', 1.5);
            end;
            hold on;
            plot(medianTable(:, 1), medianTable(:, currentIndex), 'Color', [1 0 0], 'LineWidth', 1.5);
            plot(meanTable(:, 1), meanTable(:, currentIndex), '--', 'Color', [1 0 0], 'LineWidth', 1.5);
            set(h, 'FontSize', 14, 'LineWidth', 1.5);
            xlabel(h, 'time point'); ylabel(h, [label unit]);
            currentFile = [resultsDirectory filesep '' header '.png']; saveas(h, currentFile);
            currentFrame = imread(currentFile); imwrite(currentFrame, currentFile, 'png');
            close(fh);
            % close;
        end;
end;

if averaging == 0
    lookUpTable = meanTable;
elseif averaging == 1
    lookUpTable = medianTable;
end;

if readFactors
    switch processingMode
        case 0
            lookUpTable(lookUpTable(:, 4) >= 1, 5) = 2;
            lookUpTable(lookUpTable(:, 4) < 1, 5)  = 1;
            lookUpTable(lookUpTable(:, 4) < 1, 4)  = 1 ./ lookUpTable(lookUpTable(:, 4) < 1, 4);
            
            lookUpTable(lookUpTable(:, 8) >= 1, 9) = 2;
            lookUpTable(lookUpTable(:, 8) < 1, 9)  = 1;
            lookUpTable(lookUpTable(:, 8) < 1, 8)  = 1 ./ lookUpTable(lookUpTable(:, 8) < 1, 8);
            
            lookUpTable(lookUpTable(:, 13) >= 1, 14) = 2;
            lookUpTable(lookUpTable(:, 13) < 1, 14)  = 1;
            lookUpTable(lookUpTable(:, 13) < 1, 13)  = 1 ./ lookUpTable(lookUpTable(:, 13) < 1, 13);
        case 1
            lookUpTable(lookUpTable(:, 4) >= 1, 5) = 2;
            lookUpTable(lookUpTable(:, 4) < 1, 5)  = 1;
            lookUpTable(lookUpTable(:, 4) < 1, 4)  = 1 ./ lookUpTable(lookUpTable(:, 4) < 1, 4);
        case 2
            lookUpTable(lookUpTable(:, 5) >= 1, 6) = 2;
            lookUpTable(lookUpTable(:, 5) < 1, 6)  = 1;
            lookUpTable(lookUpTable(:, 5) < 1, 5)  = 1 ./ lookUpTable(lookUpTable(:, 5) < 1, 5);
    end;
end;

save([resultsDirectory filesep 'lookUpTable.mat'], 'lookUpTable');

% warning, if data is missing
missingDataFlag = ~isempty(find(isnan(lookUpTable), 1));
if missingDataFlag
    disp(' ');
    disp('WARNING: multiFuse reference data are incomplete!');
    disp(' ');
end;