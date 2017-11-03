function collectProjections(...
    timepoints, outputString, outputID, filterID, ...
    specimen, cameras, channels, folder, ...
    inputType, outputType)

% -----------------------------------------------------------------------------------------------
% | Collection of image data projections                                                        |
% |                                                                                             |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2014                          |
% | Email: kellerp@janelia.hhmi.org                                                             |
% |                                                                                             |
% | Utilizes optimization modules and functions by Fernando Amat, HHMI/Janelia Research Campus: |
% | readKLBstack.mexw64                                                                         |
% | writeKLBstack.mexw64                                                                        |
% -----------------------------------------------------------------------------------------------

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

if ~isempty(filterID) && ~strcmp(filterID, 'corrected')
    xyBProjections = [];
    xyFProjections = [];
end;
xyProjections  = [];
xzProjections  = [];
yzProjections  = [];

if ~isempty(folder)
    mkdir(folder);
    folder = [folder filesep ''];
end;

stackInitFlag  = 0;

if length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion)');
end;

specimenHeader = ['SPM' num2str(specimen, '%.2d')];

disp(' ');

disp('---')
disp('collecting projections');
disp('---')

for currentTP = timepoints
    outputFolder = [outputString '.TM' num2str(currentTP, '%.6d') outputID];
    outputHeader = [outputFolder '/SPM' num2str(specimen, '%.2d') '_TM' num2str(currentTP, '%.6d') configurationString];
    currentIndex = find(timepoints == currentTP, 1);
    
    if ~isempty(filterID) && ~strcmp(filterID, 'corrected')
        blendedStackProjXYBName = [outputHeader '.fusedStack_xyBProjection.' filterID inputExtension];
        blendedStackProjXYFName = [outputHeader '.fusedStack_xyFProjection.' filterID inputExtension];
        blendedStackProjXYName  = [outputHeader '.fusedStack_xyProjection.' filterID inputExtension];
        blendedStackProjXZName  = [outputHeader '.fusedStack_xzProjection.' filterID inputExtension];
        blendedStackProjYZName  = [outputHeader '.fusedStack_yzProjection.' filterID inputExtension];

        if exist(blendedStackProjXYBName, 'file') == 2 && exist(blendedStackProjXYFName, 'file') == 2 && exist(blendedStackProjXYName, 'file') == 2 && exist(blendedStackProjXZName, 'file') == 2 && exist(blendedStackProjYZName, 'file') == 2
            if stackInitFlag == 0
                switch inputType
                    case 0
                        headerInformationXYB = readKLBheader(blendedStackProjXYBName);
                        imageDimensionsXYB = headerInformationXYB.xyzct(1:3);
                        
                        headerInformationXYF = readKLBheader(blendedStackProjXYFName);
                        imageDimensionsXYF = headerInformationXYF.xyzct(1:3);
                        
                        headerInformationXY = readKLBheader(blendedStackProjXYName);
                        imageDimensionsXY = headerInformationXY.xyzct(1:3);
                        
                        headerInformationXZ = readKLBheader(blendedStackProjXZName);
                        imageDimensionsXZ = headerInformationXZ.xyzct(1:3);
                        
                        headerInformationYZ = readKLBheader(blendedStackProjYZName);
                        imageDimensionsYZ = headerInformationYZ.xyzct(1:3);
                    case 1
                        [imageDimensionsXYB, bitDepth] = readJP2header(blendedStackProjXYBName);
                        
                        [imageDimensionsXYF, bitDepth] = readJP2header(blendedStackProjXYFName);
                        
                        [imageDimensionsXY, bitDepth] = readJP2header(blendedStackProjXYName);
                        
                        [imageDimensionsXZ, bitDepth] = readJP2header(blendedStackProjXZName);
                        
                        [imageDimensionsYZ, bitDepth] = readJP2header(blendedStackProjYZName);
                    case 2
                        headerInformationXYB = imfinfo(blendedStackProjXYBName);
                        imageDimensionsXYB = [headerInformationXYB(1).Height headerInformationXYB(1).Width numel(headerInformationXYB)];
                        
                        headerInformationXYF = imfinfo(blendedStackProjXYFName);
                        imageDimensionsXYF = [headerInformationXYF(1).Height headerInformationXYF(1).Width numel(headerInformationXYF)];
                        
                        headerInformationXY = imfinfo(blendedStackProjXYName);
                        imageDimensionsXY = [headerInformationXY(1).Height headerInformationXY(1).Width numel(headerInformationXY)];
                        
                        headerInformationXZ = imfinfo(blendedStackProjXZName);
                        imageDimensionsXZ = [headerInformationXZ(1).Height headerInformationXZ(1).Width numel(headerInformationXZ)];
                        
                        headerInformationYZ = imfinfo(blendedStackProjYZName);
                        imageDimensionsYZ = [headerInformationYZ(1).Height headerInformationYZ(1).Width numel(headerInformationYZ)];
                end;
                
                xyBProjections = zeros(imageDimensionsXYB(1), imageDimensionsXYB(2), numel(timepoints), 'uint16');
                xyFProjections = zeros(imageDimensionsXYF(1), imageDimensionsXYF(2), numel(timepoints), 'uint16');
                xyProjections = zeros(imageDimensionsXY(1), imageDimensionsXY(2), numel(timepoints), 'uint16');
                xzProjections = zeros(imageDimensionsXZ(1), imageDimensionsXZ(2), numel(timepoints), 'uint16');
                yzProjections = zeros(imageDimensionsYZ(1), imageDimensionsYZ(2), numel(timepoints), 'uint16');
                
                stackInitFlag = 1;
            end;
            
            xyBProjections(:, :, currentIndex) = readImage(blendedStackProjXYBName);
            xyFProjections(:, :, currentIndex) = readImage(blendedStackProjXYFName);
            xyProjections(:, :, currentIndex)  = readImage(blendedStackProjXYName);
            xzProjections(:, :, currentIndex)  = readImage(blendedStackProjXZName);
            yzProjections(:, :, currentIndex)  = readImage(blendedStackProjYZName);
        else
            disp(['missing projections: TM' num2str(currentTP, '%.6d')]);
        end;
    else
        if isempty(filterID)
            blendedStackProjXYName  = [outputHeader '.fusedStack_xyProjection' inputExtension];
            blendedStackProjXZName  = [outputHeader '.fusedStack_xzProjection' inputExtension];
            blendedStackProjYZName  = [outputHeader '.fusedStack_yzProjection' inputExtension];
        else
            blendedStackProjXYName  = [outputHeader '.fusedStack_xyProjection.' filterID inputExtension];
            blendedStackProjXZName  = [outputHeader '.fusedStack_xzProjection.' filterID inputExtension];
            blendedStackProjYZName  = [outputHeader '.fusedStack_yzProjection.' filterID inputExtension];
        end;

        if exist(blendedStackProjXYName, 'file') == 2 && exist(blendedStackProjXZName, 'file') == 2 && exist(blendedStackProjYZName, 'file') == 2
            if stackInitFlag == 0
                switch inputType
                    case 0
                        headerInformationXY = readKLBheader(blendedStackProjXYName);
                        imageDimensionsXY = headerInformationXY.xyzct(1:3);
                        
                        headerInformationXZ = readKLBheader(blendedStackProjXZName);
                        imageDimensionsXZ = headerInformationXZ.xyzct(1:3);
                        
                        headerInformationYZ = readKLBheader(blendedStackProjYZName);
                        imageDimensionsYZ = headerInformationYZ.xyzct(1:3);
                    case 1
                        [imageDimensionsXY, bitDepth] = readJP2header(blendedStackProjXYName);
                        
                        [imageDimensionsXZ, bitDepth] = readJP2header(blendedStackProjXZName);
                        
                        [imageDimensionsYZ, bitDepth] = readJP2header(blendedStackProjYZName);
                    case 2
                        headerInformationXY = imfinfo(blendedStackProjXYName);
                        imageDimensionsXY = [headerInformationXY(1).Height headerInformationXY(1).Width numel(headerInformationXY)];
                        
                        headerInformationXZ = imfinfo(blendedStackProjXZName);
                        imageDimensionsXZ = [headerInformationXZ(1).Height headerInformationXZ(1).Width numel(headerInformationXZ)];
                        
                        headerInformationYZ = imfinfo(blendedStackProjYZName);
                        imageDimensionsYZ = [headerInformationYZ(1).Height headerInformationYZ(1).Width numel(headerInformationYZ)];
                end;
                
                xyProjections = zeros(imageDimensionsXY(1), imageDimensionsXY(2), numel(timepoints), 'uint16');
                xzProjections = zeros(imageDimensionsXZ(1), imageDimensionsXZ(2), numel(timepoints), 'uint16');
                yzProjections = zeros(imageDimensionsYZ(1), imageDimensionsYZ(2), numel(timepoints), 'uint16');
                
                stackInitFlag = 1;
            end;
            
            xyProjections(:, :, currentIndex) = readImage(blendedStackProjXYName);
            xzProjections(:, :, currentIndex) = readImage(blendedStackProjXZName);
            yzProjections(:, :, currentIndex) = readImage(blendedStackProjYZName);
        else
            disp(['missing projections: TM' num2str(currentTP, '%.6d')]);
        end;
    end;
    
    if mod(find(timepoints == currentTP, 1), 50) == 0
        disp(['number of time points read from disk: ' num2str(find(timepoints == currentTP, 1))]);
    end;
end;

disp(' ');

if stackInitFlag
    disp('---')
    disp('saving projections');
    disp('---')

    if ~isempty(filterID)
        filterID = ['.' filterID];
    end;

    if length(cameras) == 2 && length(channels) == 2
        outputString = ['CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') outputID filterID]; % 4-view fusion
    elseif length(cameras) == 1 && length(channels) == 2
        outputString = ['CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d') outputID filterID]; % 2-view channel fusion
    elseif length(cameras) == 2 && length(channels) == 1
        outputString = ['CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') outputID filterID]; % 2-view camera fusion
    end;

    if ~isempty(filterID) && ~strcmp(filterID, '.corrected')
        writeImage(xyBProjections, [folder specimenHeader '_xyBProjections_' outputString outputExtension]);
        writeImage(xyFProjections, [folder specimenHeader '_xyFProjections_' outputString outputExtension]);
    end;
    writeImage(xyProjections, [folder specimenHeader '_xyProjections_' outputString outputExtension]);
    writeImage(xzProjections, [folder specimenHeader '_xzProjections_' outputString outputExtension]);
    writeImage(yzProjections, [folder specimenHeader '_yzProjections_' outputString outputExtension]);
else
     disp('no image data found at specified location');
end;

disp(' ');

end