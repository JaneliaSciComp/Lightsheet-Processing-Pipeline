folders = {...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_180.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_250.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_300.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_400.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_420.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Training.Timepoint_480.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Testing.Timepoint_120.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Testing.Timepoint_240.presumedFalsePositives\ProjectionsXY';...
    'E:\Mouse Development\Division Detection\MK5 Inspection Large Radius\Testing.Timepoint_360.presumedFalsePositives\ProjectionsXY'};

for f = 1:numel(folders)
    inputFolder = folders{f};
    outputFolder = [folders{f} '.normalized'];
    
    disp(['populating ' outputFolder]);
    
    if exist(outputFolder, 'dir') ~= 7
        mkdir(outputFolder);
    end;
    
    files = dir([inputFolder '\*.tif']);
    
    for i = 1:numel(files)
        stack = double(readImage([inputFolder '\' files(i).name]));
        stack = uint16(stack .* ((2^16 - 1) / max(stack(:))));
        for z = 1:size(stack, 3)
            if z == 1
                writeMode = 'overwrite';
            else
                writeMode = 'append';
            end;
            imwrite(stack(:, :, z), [outputFolder '\' files(i).name], 'Compression', 'none', 'WriteMode', writeMode);
        end;
    end;
end;