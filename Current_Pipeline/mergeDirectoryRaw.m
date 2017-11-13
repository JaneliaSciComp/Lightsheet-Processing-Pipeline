%% parameters

sourceDirectory  = ['D:' filesep 'Temp' filesep 'Dme_L2_57C10-GCaMP641-0_20140104_143744' filesep 'SPM00' filesep 'TM?????' filesep 'ANG000'];
sourceTimepoints = 0:10;

targetDirectory  = ['D:' filesep 'Temp' filesep 'Dme_L1_57C10-GCaMP641_0_20140104_114246' filesep 'SPM00' filesep 'TM?????' filesep 'ANG000'];
targetTimepoints = 11:21;

stampDigits      = 5;

%% main loop

disp(' ');

if numel(sourceTimepoints) ~= numel(targetTimepoints)
    error('Inconsistent number of time points');
end;

nTimepoints = numel(sourceTimepoints);
inputDatabase = cell(nTimepoints, 1);
outputDatabase = cell(nTimepoints, 1);

for n = 1:nTimepoints
    inputDatabase{n} = [sourceDirectory filesep ''];
    positions = strfind(inputDatabase{n}, repmat('?', [1 stampDigits]));
    precision = ['%.' num2str(stampDigits) 'd'];
    for k = 1:numel(positions)
        inputDatabase{n}(positions(k):(positions(k) + stampDigits - 1)) = num2str(sourceTimepoints(n), precision);
    end;
    
    outputDatabase{n} = [targetDirectory filesep ''];
    positions = strfind(outputDatabase{n}, repmat('?', [1 stampDigits]));
    precision = ['%.' num2str(stampDigits) 'd'];
    for k = 1:numel(positions)
        outputDatabase{n}(positions(k):(positions(k) + stampDigits - 1)) = num2str(targetTimepoints(n), precision);
    end;
end;

for n = 1:nTimepoints
    disp(['renaming source time point ' num2str(sourceTimepoints(n), '%.6d') ' to target time point ' num2str(targetTimepoints(n), '%.6d')]);
    
    if ~exist(outputDatabase{n}, 'dir')
        mkdir(outputDatabase{n});
    end;
    
    xmlFiles = dir([inputDatabase{n}(1:(end - 7)) '*.xml']);
    
    if ~isempty(xmlFiles)
        for x = 1:numel(xmlFiles)
            sourceFileFullPath = [inputDatabase{n}(1:(end - 7)) xmlFiles(x).name];
            targetFileFullPath = [outputDatabase{n}(1:(end - 7)) xmlFiles(x).name];
            [status, cmdout] = system(['move ' sourceFileFullPath ' ' targetFileFullPath]);
        end;
    end;
    
    sourceFiles = dir(inputDatabase{n});
    
    for i = 1:numel(sourceFiles)
        if ~strcmp(sourceFiles(i).name, '.') && ~strcmp(sourceFiles(i).name, '..')
            sourceFileName = sourceFiles(i).name;
            targetFileName = sourceFiles(i).name;
            positions = strfind(targetFileName, ['TM' num2str(sourceTimepoints(n), ['%.' num2str(stampDigits) 'd'])]);
            for k = 1:numel(positions)
                targetFileName(positions(k):(positions(k) + stampDigits + 1)) = ['TM' num2str(targetTimepoints(n), ['%.' num2str(stampDigits) 'd'])];
            end;
            
            sourceFileFullPath = [inputDatabase{n} sourceFileName];
            targetFileFullPath = [outputDatabase{n} targetFileName];
            [status, cmdout] = system(['move ' sourceFileFullPath ' ' targetFileFullPath]);
        end;
    end;
end;