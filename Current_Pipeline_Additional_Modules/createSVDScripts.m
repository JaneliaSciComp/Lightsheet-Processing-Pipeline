%% Global configuration

createShellScript = 1;
createXMLs        = 1;

includeMVDJobs    = 0;
includeSVDJobs    = 1;

checkForResults   = 1;

%% Configuration of shell script

scriptName = 'runSVD.ps1';

binaryFile = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'Binaries' filesep 'bin' filesep 'main_multiviewDeconvLR_multiGPU_blocksZ.exe';
xmlFolder  = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered' filesep 'XMLs';

%% Configuration of XML file generation

timePoints = 0:278;

stackFolder  = 'V:' filesep 'SV1' filesep 'KM_15-08-10' filesep 'Mmu_E1_mKate2_20150810_160708.corrected.registered';
stackHeaders = {...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM01_TM';...
    'SPM01_TM'};
stackFooters = {...
    '_CM00_CHN00.affine.trsf.cropped.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb';...
    '_CM00_CHN00.affine.trsf.cropped.klb';...
    '_CM01_CHN00.affine.trsf.cropped.klb'};

psfTypeArray    = {'PSF_YW_2.2_5_'; 'PSF_ZZ_2.5_8_'};
psfNameArray    = {'SmallPSF'; 'LargePSF'};
iterationArray  = [20 50];

psfFolderHeader = 'PSF';
psfHeaders      = {...
    'SPM00_TM';...
    'SPM00_TM';...
    'SPM01_TM';...
    'SPM01_TM'};
psfFooters      = {...
    '_CM00_CHN00.trimmed.klb';...
    '_CM01_CHN00.trimmed.klb';...
    '_CM00_CHN00.trimmed.klb';...
    '_CM01_CHN00.trimmed.klb'};

lambda       = 0;
background   = 100;
verbosity    = 0;
zBlockSize   = 1024;
saveUINT16   = 0;

%% Shell script generation

if createShellScript
    disp(' ');
    disp('Creating shell script');
    
    fid = fopen(scriptName, 'w');
    
    for t = 1:numel(timePoints)
        for i = 1:numel(psfTypeArray)
            outputStackName = [stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{1} num2str(timePoints(t), '%.6d') stackFooters{1} ...
                '_dec_LR_multiGPU_MVD_' psfNameArray{i} '_iter' num2str(iterationArray(i)) '_lambdaTV' num2str(lambda*10^6, '%.6d') '.klb'];
            
            if includeMVDJobs && (~checkForResults || (exist(outputStackName, 'file') ~= 2 && exist([outputStackName(1:(end-3)) 'uint16.klb'], 'file') ~= 2))
                disp(['* Including MVD XML for time point ' num2str(timePoints(t)) ' and PSF type "' psfNameArray{i} '"']);
                fwrite(fid, ['&"' binaryFile '" "' xmlFolder filesep 'TM' num2str(timePoints(t), '%.3d') '_MVD_' psfNameArray{i} '.xml"']); fprintf(fid, '\n');
            end;
            
            if includeSVDJobs
                for s = 1:numel(stackHeaders)
                    outputStackName = [stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{s} num2str(timePoints(t), '%.6d') stackFooters{s} ...
                        '_dec_LR_multiGPU_SVD_' psfNameArray{i} '_iter' num2str(iterationArray(i)) '_lambdaTV' num2str(lambda*10^6, '%.6d') '.klb'];
                    
                    if ~checkForResults || (exist(outputStackName, 'file') ~= 2 && exist([outputStackName(1:(end-3)) 'uint16.klb'], 'file') ~= 2)
                        disp(['* Including SVD XML for view ' num2str(s) ' at time point ' num2str(timePoints(t)) ' and PSF type "' psfNameArray{i} '"']);
                        fwrite(fid, ['&"' binaryFile '" "' xmlFolder filesep 'TM' num2str(timePoints(t), '%.3d') '_SVD' num2str(s) '_' psfNameArray{i} '.xml"']); fprintf(fid, '\n');
                    end;
                end;
            end;
        end;
    end;
    
    fclose(fid);
end;

%% Processing loop for XML file generation

if exist(xmlFolder, 'dir') ~= 7
    mkdir(xmlFolder);
end;

if createXMLs
    disp(' ');
    disp('Creating XML files');
    
    indentityMatrixString = ['A="'...
        '1.000000000000 0.000000000000 0.000000000000 0.000000000000 '...
        '0.000000000000 1.000000000000 0.000000000000 0.000000000000 '...
        '0.000000000000 0.000000000000 1.000000000000 0.000000000000 '...
        '0.000000000000 0.000000000000 0.000000000000 1.000000000000"'];
    
    for t = 1:numel(timePoints)
        for i = 1:numel(psfTypeArray)
            outputStackName = [stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{1} num2str(timePoints(t), '%.6d') stackFooters{1} ...
                '_dec_LR_multiGPU_MVD_' psfNameArray{i} '_iter' num2str(iterationArray(i)) '_lambdaTV' num2str(lambda*10^6, '%.6d') '.klb'];
            
            if includeMVDJobs && (~checkForResults || (exist(outputStackName, 'file') ~= 2 && exist([outputStackName(1:(end-3)) 'uint16.klb'], 'file') ~= 2))
                disp(['* Creating MVD XML for time point ' num2str(timePoints(t)) ' and PSF type "' psfNameArray{i} '"']);
                
                outputName = [xmlFolder filesep 'TM' num2str(timePoints(t), '%.3d') '_MVD_' psfNameArray{i} '.xml'];
                
                fid = fopen(outputName, 'w');
                
                fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>'); fprintf(fid, '\n');
                fwrite(fid, '<document>'); fprintf(fid, '\n');
                
                for s = 1:numel(stackHeaders)
                    fwrite(fid, [...
                        '<view imgFilename="' stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{s} num2str(timePoints(t), '%.6d') stackFooters{s}  ...
                        '" psfFilename="' stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' psfFolderHeader filesep '' psfTypeArray{i} psfHeaders{s} num2str(timePoints(t), '%.6d') psfFooters{s} ...
                        '" ' indentityMatrixString '></view>']); fprintf(fid, '\n');
                end;
                
                fwrite(fid, ['<deconvolution lambdaTV="' num2str(lambda, '%.6f') ...
                    '" numIter="' num2str(iterationArray(i)) ...
                    '" imBackground="' num2str(background, '%.6f') ...
                    '" verbose="' num2str(verbosity) ...
                    '" blockZsize="' num2str(zBlockSize) ...
                    '" prefix="MVD_' psfNameArray{i} ...
                    '" SaveAsUINT16="' num2str(saveUINT16) '"/>']); fprintf(fid, '\n');
                
                fwrite(fid, '</document>'); fprintf(fid, '\n');
                
                fclose(fid);
            end;
            
            if includeSVDJobs
                for s = 1:numel(stackHeaders)
                    outputStackName = [stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{s} num2str(timePoints(t), '%.6d') stackFooters{s} ...
                        '_dec_LR_multiGPU_SVD_' psfNameArray{i} '_iter' num2str(iterationArray(i)) '_lambdaTV' num2str(lambda*10^6, '%.6d') '.klb'];
                    
                    if ~checkForResults || (exist(outputStackName, 'file') ~= 2 && exist([outputStackName(1:(end-3)) 'uint16.klb'], 'file') ~= 2)
                        disp(['* Creating SVD XML for view ' num2str(s) ' at time point ' num2str(timePoints(t)) ' and PSF type "' psfNameArray{i} '"']);
                        
                        outputName = [xmlFolder filesep 'TM' num2str(timePoints(t), '%.3d') '_SVD' num2str(s) '_' psfNameArray{i} '.xml'];
                        
                        fid = fopen(outputName, 'w');
                        
                        fwrite(fid, '<?xml version="1.0" encoding="utf-8"?>'); fprintf(fid, '\n');
                        fwrite(fid, '<document>'); fprintf(fid, '\n');
                        
                        fwrite(fid, [...
                            '<view imgFilename="' stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' stackHeaders{s} num2str(timePoints(t), '%.6d') stackFooters{s}  ...
                            '" psfFilename="' stackFolder filesep 'TM' num2str(timePoints(t), '%.6d') filesep '' psfFolderHeader filesep '' psfTypeArray{i} psfHeaders{s} num2str(timePoints(t), '%.6d') psfFooters{s} ...
                            '" ' indentityMatrixString '></view>']); fprintf(fid, '\n');
                        
                        fwrite(fid, ['<deconvolution lambdaTV="' num2str(lambda, '%.6f') ...
                            '" numIter="' num2str(iterationArray(i)) ...
                            '" imBackground="' num2str(background, '%.6f') ...
                            '" verbose="' num2str(verbosity) ...
                            '" blockZsize="' num2str(zBlockSize) ...
                            '" prefix="SVD_' psfNameArray{i} ...
                            '" SaveAsUINT16="' num2str(saveUINT16) '"/>']); fprintf(fid, '\n');
                        
                        fwrite(fid, '</document>'); fprintf(fid, '\n');
                        
                        fclose(fid);
                    end;
                end;
            end;
        end;
    end;
end;

disp(' ');