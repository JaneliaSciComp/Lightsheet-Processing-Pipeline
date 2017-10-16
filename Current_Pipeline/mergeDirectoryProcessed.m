%% parameters

sourceDirectory  = 'D:\Temp\Dme_L2_57C10-GCaMP641-0_20140104_143744.corrected\SPM00\TM??????';
sourceTimepoints = 0:10;

targetDirectory  = 'D:\Temp\Dme_L1_57C10-GCaMP641_0_20140104_114246.corrected\SPM00\TM??????';
targetTimepoints = 11:21;

stampDigits      = 6;

%% main loop

disp(' ');

if numel(sourceTimepoints) ~= numel(targetTimepoints)
    error('Inconsistent number of time points');
end;

nTimepoints = numel(sourceTimepoints);
inputDatabase = cell(nTimepoints, 1);
outputDatabase = cell(nTimepoints, 1);

for n = 1:nTimepoints
    inputDatabase{n} = [sourceDirectory '\'];
    positions = strfind(inputDatabase{n}, repmat('?', [1 stampDigits]));
    precision = ['%.' num2str(stampDigits) 'd'];
    for k = 1:numel(positions)
        inputDatabase{n}(positions(k):(positions(k) + stampDigits - 1)) = num2str(sourceTimepoints(n), precision);
    end;
    
    outputDatabase{n} = [targetDirectory '\'];
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
    
    projectionFiles = dir([inputDatabase{n}(1:(end - 16)) '.projections\*TM' num2str(sourceTimepoints(n), '%.6d') '*.*']);
    if ~isempty(projectionFiles)
        for i = 1:numel(projectionFiles)
            sourceFileName = projectionFiles(i).name;
            targetFileName = projectionFiles(i).name;
            positions = strfind(targetFileName, ['TM' num2str(sourceTimepoints(n), ['%.' num2str(stampDigits) 'd'])]);
            for k = 1:numel(positions)
                targetFileName(positions(k):(positions(k) + stampDigits + 1)) = ['TM' num2str(targetTimepoints(n), ['%.' num2str(stampDigits) 'd'])];
            end;
            
            sourceFileFullPath = [inputDatabase{n}(1:(end - 16)) '.projections\' sourceFileName];
            targetFileFullPath = [outputDatabase{n}(1:(end - 16)) '.projections\' targetFileName];
            [status, cmdout] = system(['move ' sourceFileFullPath ' ' targetFileFullPath]);
        end;
    end;
end;