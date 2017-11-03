dataRoot   = 'W:' filesep 'SiMView1' filesep '13-12-16' filesep 'Mmu_E1_CAGTAG1_01_20131216_164117.corrected' filesep 'SPM00';
identifier = 'Mmu_E1_CAGTAG1_01_20131216_164117';
specimen   = 0;
timePoints = 0:254;
cameras    = 0:1;
channels   = 0:1;

inputType  = 0;     % 0: input data in KLB format
                    % 1: input data in JP2 format
                    % 2: input data in TIF format

missingStacksTable = zeros(numel(timePoints), numel(cameras) * numel(channels) + 1);
missingStacksTable(:, 1) = timePoints;

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

for t = timePoints
    for c = cameras
        for h = channels
            stackFile = [dataRoot filesep 'TM' num2str(t, '%.6d') filesep 'SPM' num2str(specimen, '%.2d') ...
                '_TM' num2str(t, '%.6d') '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') inputExtension];
            xyMaskFile = [stackFile(1:(end - 3)) 'xyMask' inputExtension];
            xzMaskFile = [stackFile(1:(end - 3)) 'xzMask' inputExtension];
            
            stackFlag = 1;
            xyMaskFlag = 1;
            xzMaskFlag = 1;
            
            if ~exist(stackFile, 'file')
                disp(['missing stack file: ' stackFile]);
                stackFlag = 0;
                missingStacksTable(find(timePoints == t, 1), (find(cameras == c, 1) - 1) * 2 + find(channels == h, 1) + 1) = 1;
            end;
            if ~exist(xyMaskFile, 'file')
                disp(['missing xyMask file: ' stackFile]);
                xyMaskFlag = 0;
            end;
            if ~exist(xzMaskFile, 'file')
                disp(['missing xzMask file: ' stackFile]);
                xzMaskFlag = 0;
            end;
            
            if stackFlag && (~xyMaskFlag || ~xzMaskFlag)
                stack = readImage(stackFile) > 0;
                
                if ~xzMaskFlag
                    xzSliceMask = zeros(size(stack, 1), size(stack, 3), 'uint16');
                    coordinateMask = repmat(1:size(stack, 2), [size(stack, 1) 1 size(stack, 3)]);
                    coordinateMask(stack == 0) = NaN;
                    xzSliceMask = uint16(round(squeeze(nanmean(coordinateMask, 2))));
                    clear coordinateMask;
                    
                    writeImage(xzSliceMask, xzMaskFile);
                end;
                
                if ~xyMaskFlag
                    xySliceMask = zeros(size(stack, 1), size(stack, 2), 'uint16');
                    coordinateMask = repmat(reshape(1:size(stack, 3), [1 1 size(stack, 3)]), [size(stack, 1) size(stack, 2) 1]);
                    coordinateMask(stack == 0) = NaN;
                    xySliceMask = uint16(round(squeeze(nanmean(coordinateMask, 3))));
                    clear coordinateMask;
                    
                    writeImage(xySliceMask, xyMaskFile);
                end;
            end;
        end;
    end;
end;

save([identifier '.missingStacksTable.mat'], 'missingStacksTable');

missingStacks = sum(missingStacksTable(:, 2:(numel(cameras) * numel(channels) + 1)), 2);
missingStacks = missingStacksTable(missingStacks > 0, 1);

save([identifier '.missingStacks.txt'], 'missingStacks', '-ascii');