validDimensions = [128 243 256 512 729 1024 2048 2187 4096 6561 8192 16384];
timepointsToExclude = [];

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        else
            allDoneFlag = ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb'], 'file') == 2;
            if allDoneFlag
                candidates(i) = [];
            end;
        end;
    end;
end;

for i = 1:numel(candidates)
    timepoint = str2num(candidates(i).name(3:end));
    
    if isempty(find(timepointsToExclude == timepoint, 1))
        disp(['padding time point ' num2str(timepoint)]);
        
        stack = readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb']);
        [xSize, ySize, zSize] = size(stack);
        xSizePadded = validDimensions(find(validDimensions >= xSize, 1, 'first'));
        ySizePadded = validDimensions(find(validDimensions >= ySize, 1, 'first'));
        zSizePadded = validDimensions(find(validDimensions >= zSize, 1, 'first'));
        
        if xSize == xSizePadded && ySize == ySizePadded && zSize == zSizePadded
            copyfile(...
                [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'],...
                [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb']);
            copyfile(...
                [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'],...
                [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb']);
            copyfile(...
                [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb'],...
                [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb']);
            copyfile(...
                [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb'],...
                [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb']);
        else
            background = min(stack(:));
            if background ~= 100
                disp(['Warning: background level of ' num2str(background) ' detected for time point ' num2str(timepoint)]);
            end;
            
            paddedStack = ones(xSizePadded, ySizePadded, zSizePadded, 'uint16');
            paddedStack(...
                (round((xSizePadded - xSize) / 2) + 1):(round((xSizePadded - xSize) / 2) + xSize),...
                (round((ySizePadded - ySize) / 2) + 1):(round((ySizePadded - ySize) / 2) + ySize),...
                (round((zSizePadded - zSize) / 2) + 1):(round((zSizePadded - zSize) / 2) + zSize)) = stack;
            writeImage(paddedStack, [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb']);
            
            stack = readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb']);
            paddedStack = ones(xSizePadded, ySizePadded, zSizePadded, 'uint16');
            paddedStack(...
                (round((xSizePadded - xSize) / 2) + 1):(round((xSizePadded - xSize) / 2) + xSize),...
                (round((ySizePadded - ySize) / 2) + 1):(round((ySizePadded - ySize) / 2) + ySize),...
                (round((zSizePadded - zSize) / 2) + 1):(round((zSizePadded - zSize) / 2) + zSize)) = stack;
            writeImage(paddedStack, [candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb']);
            
            stack = readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.klb']);
            paddedStack = ones(xSizePadded, ySizePadded, zSizePadded, 'uint16');
            paddedStack(...
                (round((xSizePadded - xSize) / 2) + 1):(round((xSizePadded - xSize) / 2) + xSize),...
                (round((ySizePadded - ySize) / 2) + 1):(round((ySizePadded - ySize) / 2) + ySize),...
                (round((zSizePadded - zSize) / 2) + 1):(round((zSizePadded - zSize) / 2) + zSize)) = stack;
            writeImage(paddedStack, [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf.cropped.padded.klb']);
            
            stack = readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.klb']);
            paddedStack = ones(xSizePadded, ySizePadded, zSizePadded, 'uint16');
            paddedStack(...
                (round((xSizePadded - xSize) / 2) + 1):(round((xSizePadded - xSize) / 2) + xSize),...
                (round((ySizePadded - ySize) / 2) + 1):(round((ySizePadded - ySize) / 2) + ySize),...
                (round((zSizePadded - zSize) / 2) + 1):(round((zSizePadded - zSize) / 2) + zSize)) = stack;
            writeImage(paddedStack, [candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf.cropped.padded.klb']);
        end;
    else
        disp(['skipping time point ' num2str(timepoint)]);
    end;
end;