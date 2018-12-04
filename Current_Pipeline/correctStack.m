function correctStack(inputDatabase, outputDatabase, dataType, timepoints, dimensions, dimensionsMax, dimensionsDeltas, ...
    correctIntensity, intensityBackgrounds, intensityFactors, correctDrift, inputType, outputType, ...
    percentile, referenceROI, xOffsets, yOffsets, zOffsets, jobMemory, t, memoryEstimate)

timepoint = timepoints(t);

version = 1.01;

configuration = cell(21, 1);

configuration{1}  = version;          configuration{2}  = timepoint;            configuration{3}  = inputDatabase;    configuration{4}  = outputDatabase;
configuration{5}  = dataType;         configuration{6}  = timepoints;
configuration{7}  = dimensions;       configuration{8}  = dimensionsMax;        configuration{9}  = dimensionsDeltas;
configuration{10} = correctIntensity; configuration{11} = intensityBackgrounds; configuration{12} = intensityFactors;
configuration{13} = correctDrift;     configuration{14} = inputType;            configuration{15} = outputType;
configuration{16} = percentile;       configuration{17} = referenceROI;
configuration{18} = xOffsets;         configuration{19} = yOffsets;             configuration{20} = zOffsets;
configuration{21} = [jobMemory(1) memoryEstimate];

if ~exist(outputDatabase{t, 1}, 'dir')
    mkdir(outputDatabase{t, 1});
end;

save([outputDatabase{t, 2} '.configuration.mat'], 'configuration');

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

currentIndex = find(dimensions(:, 1) == timepoint, 1);
mismatchFlag = ...
    (dimensionsMax(1) ~= dimensions(currentIndex, 2)) || ...
    (dimensionsMax(2) ~= dimensions(currentIndex, 3)) || ...
    (dimensionsMax(3) ~= dimensions(currentIndex, 4));

if mismatchFlag
    if dataType == 0
        if correctIntensity
            unpaddedImageXY = (readImage([inputDatabase{t} '_xyProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
            unpaddedImageXZ = (readImage([inputDatabase{t} '_xzProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
            unpaddedImageYZ = (readImage([inputDatabase{t} '_yzProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
        else
            unpaddedImageXY = readImage([inputDatabase{t} '_xyProjection' inputExtension]);
            unpaddedImageXZ = readImage([inputDatabase{t} '_xzProjection' inputExtension]);
            unpaddedImageYZ = readImage([inputDatabase{t} '_yzProjection' inputExtension]);
        end;
        
        if percentile(1) == 0
            unpaddedImageXYMin = min(unpaddedImageXY(:));
            unpaddedImageXZMin = min(unpaddedImageXZ(:));
            unpaddedImageYZMin = min(unpaddedImageYZ(:));
        else
            unpaddedImageXYMin = prctile(unpaddedImageXY(:), percentile(1));
            unpaddedImageXZMin = prctile(unpaddedImageXZ(:), percentile(1));
            unpaddedImageYZMin = prctile(unpaddedImageYZ(:), percentile(1));
        end;
        
        currentImageXY = ones(dimensionsMax(1), dimensionsMax(2), 'uint16') .* unpaddedImageXYMin;
        currentImageXZ = ones(dimensionsMax(1), dimensionsMax(3), 'uint16') .* unpaddedImageXZMin;
        currentImageYZ = ones(dimensionsMax(2), dimensionsMax(3), 'uint16') .* unpaddedImageYZMin;
        
        currentImageXY(...
            (dimensionsDeltas(currentIndex, 2) + 1):(dimensionsDeltas(currentIndex, 2) + dimensions(currentIndex, 2)), ...
            (dimensionsDeltas(currentIndex, 3) + 1):(dimensionsDeltas(currentIndex, 3) + dimensions(currentIndex, 3))) = unpaddedImageXY;
        currentImageXZ(...
            (dimensionsDeltas(currentIndex, 2) + 1):(dimensionsDeltas(currentIndex, 2) + dimensions(currentIndex, 2)), ...
            (dimensionsDeltas(currentIndex, 4) + 1):(dimensionsDeltas(currentIndex, 4) + dimensions(currentIndex, 4))) = unpaddedImageXZ;
        currentImageYZ(...
            (dimensionsDeltas(currentIndex, 3) + 1):(dimensionsDeltas(currentIndex, 3) + dimensions(currentIndex, 3)), ...
            (dimensionsDeltas(currentIndex, 4) + 1):(dimensionsDeltas(currentIndex, 4) + dimensions(currentIndex, 4))) = unpaddedImageYZ;
        
        unpaddedImageXY = [];
        unpaddedImageXZ = [];
        unpaddedImageYZ = [];
    else
        if correctIntensity
            unpaddedStack = (readImage([inputDatabase{t} inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
        else
            unpaddedStack = readImage([inputDatabase{t} inputExtension]);
        end;
        
        if percentile(1) == 0
            currentStackMin = min(unpaddedStack(1:percentile(2):end));
        else
            currentStackMin = prctile(unpaddedStack(1:percentile(2):end), percentile(1));
        end;
        
        currentStack = ones(dimensionsMax(1), dimensionsMax(2), dimensionsMax(3), 'uint16') .* currentStackMin;
        
        currentStack(...
            (dimensionsDeltas(currentIndex, 2) + 1):(dimensionsDeltas(currentIndex, 2) + dimensions(currentIndex, 2)), ...
            (dimensionsDeltas(currentIndex, 3) + 1):(dimensionsDeltas(currentIndex, 3) + dimensions(currentIndex, 3)), ...
            (dimensionsDeltas(currentIndex, 4) + 1):(dimensionsDeltas(currentIndex, 4) + dimensions(currentIndex, 4))) = unpaddedStack;
        
        unpaddedStack = [];
    end;
else
    if dataType == 0
        if correctIntensity
            currentImageXY = (readImage([inputDatabase{t} '_xyProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
            currentImageXZ = (readImage([inputDatabase{t} '_xzProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
            currentImageYZ = (readImage([inputDatabase{t} '_yzProjection' inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
        else
            currentImageXY = readImage([inputDatabase{t} '_xyProjection' inputExtension]);
            currentImageXZ = readImage([inputDatabase{t} '_xzProjection' inputExtension]);
            currentImageYZ = readImage([inputDatabase{t} '_yzProjection' inputExtension]);
        end;
    else
        if correctIntensity
            currentStack = (readImage([inputDatabase{t} inputExtension]) - intensityBackgrounds(currentIndex)) .* intensityFactors(currentIndex);
        else
            currentStack = readImage([inputDatabase{t} inputExtension]);
        end;
    end;
end;

if correctDrift
    if dataType == 0
        if ~(xOffsets(currentIndex) == 0 && yOffsets(currentIndex) == 0 && zOffsets(currentIndex) == 0)
            currentImageXY = translateImage(repmat(currentImageXY, [1 1 2]), xOffsets(currentIndex), yOffsets(currentIndex), 0);
            currentImageXY = currentImageXY(:, :, 1);
            currentImageXZ = translateImage(reshape(currentImageXZ, [size(currentImageXZ, 1) 1 size(currentImageXZ, 2)]), xOffsets(currentIndex), 0, zOffsets(currentIndex));
            currentImageXZ = reshape(currentImageXZ, [size(currentImageXZ, 1) size(currentImageXZ, 3) 1]);
            currentImageYZ = translateImage(reshape(currentImageYZ, [size(currentImageYZ, 1) 1 size(currentImageYZ, 2)]), yOffsets(currentIndex), 0, zOffsets(currentIndex));
            currentImageYZ = reshape(currentImageYZ, [size(currentImageYZ, 1) size(currentImageYZ, 3) 1]);
        end;
        if ~isempty(referenceROI)
            currentImageXY = currentImageXY(referenceROI(1, 1):referenceROI(1, 2), referenceROI(2, 1):referenceROI(2, 2));
            currentImageXZ = currentImageXZ(referenceROI(1, 1):referenceROI(1, 2), referenceROI(3, 1):referenceROI(3, 2));
            currentImageYZ = currentImageYZ(referenceROI(2, 1):referenceROI(2, 2), referenceROI(3, 1):referenceROI(3, 2));
        end;
    else
        if ~(xOffsets(currentIndex) == 0 && yOffsets(currentIndex) == 0 && zOffsets(currentIndex) == 0)
            currentStack = translateImage(currentStack, xOffsets(currentIndex), yOffsets(currentIndex), zOffsets(currentIndex));
        end;
        if ~isempty(referenceROI)
            currentStack = currentStack(referenceROI(1, 1):referenceROI(1, 2), referenceROI(2, 1):referenceROI(2, 2), referenceROI(3, 1):referenceROI(3, 2));
        end;
    end;
end;

if dataType == 0
    writeImage(currentImageXY, [outputDatabase{t, 2} '_xyProjection.corrected' outputExtension]);
    writeImage(currentImageXZ, [outputDatabase{t, 2} '_xzProjection.corrected' outputExtension]);
    writeImage(currentImageYZ, [outputDatabase{t, 2} '_yzProjection.corrected' outputExtension]);
else
    writeImage(currentStack, [outputDatabase{t, 2} '.corrected' outputExtension]);
    
    writeImage(max(currentStack, [], 3), [outputDatabase{t, 2} '_xyProjection.corrected' outputExtension]);
    writeImage(squeeze(max(currentStack, [], 2)), [outputDatabase{t, 2} '_xzProjection.corrected' outputExtension]);
    writeImage(squeeze(max(currentStack, [], 1)), [outputDatabase{t, 2} '_yzProjection.corrected' outputExtension]);
end;

end