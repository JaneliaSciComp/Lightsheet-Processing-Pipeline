%% parameters

executionMode  = 1;              % 1: unmask all KLB files detected in current folder
                                 % 2: unmask all KLB files in experiment directory structure specified by the following formatting parameters

% formatting parameters for executionMode == 2
rootFolder     = 'S:' filesep 'SiMView1' filesep '15-04-17' filesep 'Mmu_E1_mKate2_20150417_153156.corrected';
specimens      = [0 1];
timepoints     = 0:1192;
cameras        = [0 1];
channels       = 0;

% global parameters                                 
correctionType = 1;              % 1: adaptive unmasking (smooth transition to background intensity)
                                 % 2: unmasking using fixed background intensity
correctionSlot = [2 100];        % first value set to 1 or 2: first value identifies minIntensity slot used as background intensity, second value ignored
                                 % first value set to 3: second value defines background intenstiy (override minIntensity data)

fastInitialize = 1;              % (correctionType == 1 only) flag for using only surface foreground voxels when executing scatteredInterpolant
erosionRadius  = 2;              % (correctionType == 1 only) radius of erosion kernel used to find foreground voxels relevant for scatteredInterpolant
interpRadius   = 5;              % (correctionType == 1 only) size of smooth transition zone bridging foreground and background
scalingFactor  = 2.031/(6.5/16); % (correctionType == 1 only) ratio of axial-vs-lateral voxel size
downsampling   = 3;              % (correctionType == 1 only) if >1, this is the downsampling factor applied to image stack for constructing interpolation function
saveMasks      = 0;              % (correctionType == 1 only) flag for saving intermediate binary and dilated masks
saveMetadata   = 0;              % (correctionType == 1 only) flag for saving voxel count and timing information

useProfiling   = 0;              % activate memory and computation time profiling
verbose        = 0;              % flag for displaying information on current processing progress
nParallelJobs  = 4;              % number of jobs that are processed in parallel (set to 1 to disable parallel processing)

%% main loop

if useProfiling
    profile('-memory', 'on');
    setpref('profiler', 'showJitLines', 1);
    
    profile on;
end;

if executionMode == 1
    stackNames = dir('*.klb');
    inputFolder = pwd;
    outputFolder = [inputFolder filesep 'Unmasked'];
    
    if ~isempty(stackNames)
        disp(' ');
        
        if nParallelJobs == 1
            for s = 1:length(stackNames)
                maskedStackName = stackNames(s).name;
                
                unmaskStack(...
                    inputFolder, maskedStackName,...
                    correctionType, correctionSlot,...
                    fastInitialize, erosionRadius, interpRadius, scalingFactor, downsampling,...
                    saveMasks, saveMetadata,...
                    outputFolder, verbose, nParallelJobs)
                
                if verbose
                    disp(' ');
                end;
            end;
            
            if ~verbose
                disp(' ');
            end;
        else
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            matlabpool(nParallelJobs);
            disp(' ');
            
            parfor s = 1:length(stackNames)
                maskedStackName = stackNames(s).name;
                
                unmaskStack(...
                    inputFolder, maskedStackName,...
                    correctionType, correctionSlot,...
                    fastInitialize, erosionRadius, interpRadius, scalingFactor, downsampling,...
                    saveMasks, saveMetadata,...
                    outputFolder, verbose, nParallelJobs)
            end;
            
            disp(' ');
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            disp(' ');
        end;
    else
        disp(' ');
        disp('No KLB files found');
        disp(' ');
    end;
else
    disp(' ');
    missingData = [];
    
    if nParallelJobs == 1
        for t = timepoints
            for s = specimens
                for c = cameras
                    for h = channels
                        inputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                        outputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                        maskedStackName = ['SPM' num2str(s, '%.2d') '_TM' num2str(t, '%.6d') '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.klb'];
                        
                        if exist([inputFolder filesep '' maskedStackName], 'file') == 2
                            unmaskStack(...
                                inputFolder, maskedStackName,...
                                correctionType, correctionSlot,...
                                fastInitialize, erosionRadius, interpRadius, scalingFactor, downsampling,...
                                saveMasks, saveMetadata,...
                                outputFolder, verbose, nParallelJobs)
                            
                            if verbose
                                disp(' ');
                            end;
                        else
                            missingData = cat(1, missingData, [t, s, c, h]);
                            disp(['Could not find file ' maskedStackName]);
                            
                            if verbose
                                disp(' ');
                            end;
                        end;
                    end;
                end;
            end;
        end;
    else
        disp('Populating configuration table');
        configurationTable = [];
        for t = timepoints
            for s = specimens
                for c = cameras
                    for h = channels
                        inputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                        outputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                        maskedStackName = ['SPM' num2str(s, '%.2d') '_TM' num2str(t, '%.6d') '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.klb'];
                        
                        if exist([inputFolder filesep '' maskedStackName], 'file') == 2
                            configurationTable = cat(1, configurationTable, [t, s, c, h]);
                        else
                            missingData = cat(1, missingData, [t, s, c, h]);
                            disp(['Could not find file ' maskedStackName]);
                        end;
                    end;
                end;
            end;
        end;
        
        if ~isempty(configurationTable)
            disp(' ');
            if matlabpool('size') > 0
                matlabpool('close');
            end;
            matlabpool(nParallelJobs);
            disp(' ');
            
            parfor n = 1:size(configurationTable, 1)
                t = configurationTable(n, 1);
                s = configurationTable(n, 2);
                c = configurationTable(n, 3);
                h = configurationTable(n, 4);
                inputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                outputFolder = [rootFolder filesep 'SPM' num2str(s, '%.2d') filesep 'TM' num2str(t, '%.6d')];
                maskedStackName = ['SPM' num2str(s, '%.2d') '_TM' num2str(t, '%.6d') '_CM' num2str(c, '%.2d') '_CHN' num2str(h, '%.2d') '.klb'];
                
                unmaskStack(...
                    inputFolder, maskedStackName,...
                    correctionType, correctionSlot,...
                    fastInitialize, erosionRadius, interpRadius, scalingFactor, downsampling,...
                    saveMasks, saveMetadata,...
                    outputFolder, verbose, nParallelJobs)
            end;
            
            disp(' ');
            if matlabpool('size') > 0
                matlabpool('close');
            end;
        else
            disp(' ');
            disp('Configuration table is empty - aborting parallel job submission');
        end;
    end;
    
    if ~isempty(missingData)
        currentTime = clock;
        timeString = [num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
            '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
        missingDataName = [pwd filesep 'missingData.unmaskStack.' timeString '.mat'];
        save(missingDataName, 'missingData');
        disp(' ');
        disp('WARNING: not all files were found!');
        disp(['Missing files are listed in missingData.unmaskStack.' timeString '.mat']);
        disp(' ');
    else
        if ~verbose
            disp(' ');
        end;
    end;
end;

if useProfiling
    profile off;
    
    currentTime = clock;
    timeString = [num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
        '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
    if fastInitialize
        profileName = ['profiling_r' num2str(interpRadius) '_d' num2str(downsampling) '_e' num2str(erosionRadius) '.' timeString];
    else
        profileName = ['profiling_r' num2str(interpRadius) '_d' num2str(downsampling) '.' timeString];
    end;
    profsave(profile('info'), profileName);
end;