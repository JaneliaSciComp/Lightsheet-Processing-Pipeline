folders = {...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_180.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_250.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_300.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_400.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_420.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Training.Timepoint_480.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Testing.Timepoint_120.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Testing.Timepoint_240.presumedFalsePositives' filesep 'ProjectionsXY';...
    'E:' filesep 'Mouse Development' filesep 'Division Detection' filesep 'MK5 Inspection Large Radius' filesep 'Testing.Timepoint_360.presumedFalsePositives' filesep 'ProjectionsXY'};

for f = 1:numel(folders)
    inputFolder = folders{f};
    outputFolder = [folders{f} '.normalized'];
    
    disp(['populating ' outputFolder]);
    
    if exist(outputFolder, 'dir') ~= 7
        mkdir(outputFolder);
    end;
    
    files = dir([inputFolder filesep '*.tif']);
    
    for i = 1:numel(files)
        stack = double(readImage([inputFolder filesep '' files(i).name]));
        stack = uint16(stack .* ((2^16 - 1) / max(stack(:))));
        for z = 1:size(stack, 3)
            if z == 1
                writeMode = 'overwrite';
            else
                writeMode = 'append';
            end;
            imwrite(stack(:, :, z), [outputFolder filesep '' files(i).name], 'Compression', 'none', 'WriteMode', writeMode);
        end;
    end;
end;