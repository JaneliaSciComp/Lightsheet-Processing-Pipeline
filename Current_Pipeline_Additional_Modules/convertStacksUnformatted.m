inputString1  = 'V:\SV1\KM_15-08-10\ngc\WeightedXRestrictedRawCentroidTiltedZones\SPM00_TM';
inputString3  = '_CM00_CHN';
inputString4  = '.weighted_xRestrictedRawCentroidsTiltedZones.klb_dec_LR_multiGPU_FWHMpsfZ_05_iter20_lambdaTV000000.klb';
outputFolder  = 'V:\SV1\KM_15-08-10\ngc\Downsampled';
timepoints    = 0:278;
channels      = 0;
resolution    = [6.5/16 6.5/16 6.5/16];
padding       = 1;   % match z-size of all stacks to global maximum (computation of maximum depends on parameter fullSearch)
fullSearch    = 1;   % 0: assume that z-size is largest at last time point, 1: determine global maximum by analyzing all time points
downsampling  = 4;   % downsampling factor (crop to integer-multiple of this number first)
isotropicFlag = 1;   % 0: downsampling only in x and y, 1: downsampling in x, y and z
ouputFormat   = 'inr';

machineFormat = 'l'; % only required for 'inr' file format
                     % use 'l' for correct byte order in block-matching software (default)
                     % use 'b' for correct byte order in Fiji
parallelJobs  = 12;

%% main loop

resolution([1 2]) = resolution([1 2]) .* downsampling;

if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end;

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(parallelJobs);

if padding
    timeString = num2str(timepoints(end), '%.6d');
    channelString = num2str(channels(end), '%.2d');
    stackName = [inputString1 timeString inputString3 channelString inputString4];
    stackInformation = readKLBheader(stackName);
    globalZSize = stackInformation(1).xyzct(3);
    
    if fullSearch
        for t = 1:numel(timepoints)
            timeString = num2str(timepoints(t), '%.6d');
            stackName = [inputString1 timeString inputString3 channelString inputString4];
            stackInformation = readKLBheader(stackName);
            globalZSize = max(globalZSize, stackInformation(1).xyzct(3));
        end;
    end;
end;

parfor t = 1:numel(timepoints)
    for c = 1:numel(channels)
        disp(['Converting image stack at TM' num2str(timepoints(t)) ', channel ' num2str(channels(c))]);
        
        timeString = num2str(timepoints(t), '%.6d');
        channelString = num2str(channels(c), '%.2d');
        stackName = [inputString1 timeString inputString3 channelString inputString4];
        stack = readImage(stackName);
        
        [xSize, ySize, zSize] = size(stack);
        newXSize = (xSize - mod(xSize, downsampling)) / downsampling;
        newYSize = (ySize - mod(ySize, downsampling)) / downsampling;
        if padding
            newZSize = globalZSize;
        else
            newZSize = zSize;
        end;
        currentZOffset = round((newZSize - zSize) / 2);
        convertedStack = zeros(newXSize, newYSize, newZSize, 'uint16');
        
        for x = 1:newXSize
            compressedBlock = mean(stack((1:downsampling) + (x-1) * downsampling, :, :), 1);
            for y = 1:newYSize
                convertedStack(x, y, (currentZOffset + 1):(currentZOffset + zSize)) = mean(compressedBlock(1, (1:downsampling) + (y-1) * downsampling, :), 2);
            end;
        end;
        
        if isotropicFlag
            newZSize = (newZSize - mod(newZSize, downsampling)) / downsampling;
            finalStack = zeros(newXSize, newYSize, newZSize, 'uint16');
            
            for z = 1:newZSize
                finalStack(:, :, z) = mean(convertedStack(:, :, (1:downsampling) + (z-1) * downsampling), 3);
            end;
            
            convertedStack = finalStack;
        end;
        
        outputName = [outputFolder '\TM' num2str(timepoints(t), '%.3d') '_CHN' num2str(channels(c)) '.Downsampled' num2str(downsampling) '.' ouputFormat];
        
        if strcmp(ouputFormat, 'inr')
            fid = fopen(outputName, 'w');
            headerLength = 0;
            
            headerA = '#INRIMAGE-4#{';
            fwrite(fid, headerA);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerA) + 1;
            
            headerB = ['XDIM=' num2str(newXSize)];
            fwrite(fid, headerB);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerB) + 1;
            
            headerC = ['YDIM=' num2str(newYSize)];
            fwrite(fid, headerC);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerC) + 1;
            
            headerD = ['ZDIM=' num2str(newZSize)];
            fwrite(fid, headerD);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerD) + 1;
            
            headerE = 'TYPE=unsigned fixed'; % can also be float
            fwrite(fid, headerE);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerE) + 1;
            
            headerF = 'PIXSIZE=16 bits'; % can be 8, 16, 32, 64 bit
            fwrite(fid, headerF);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerF) + 1;
            
            headerG = 'CPU=decm';
            fwrite(fid, headerG);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerG) + 1;
            
            headerH = ['VX=' num2str(resolution(1))];
            fwrite(fid, headerH);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerH) + 1;
            
            headerI = ['VY=' num2str(resolution(2))];
            fwrite(fid, headerI);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerI) + 1;
            
            headerJ = ['VZ=' num2str(resolution(3))];
            fwrite(fid, headerJ);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerJ) + 1;
            
            headerK = '#GEOMETRY=CARTESIAN';
            fwrite(fid, headerK);
            fprintf(fid, '\n');
            headerLength = headerLength + length(headerK) + 1;
            
            if mod(headerLength + 4, 256) ~= 0 % total header length has to be multiple of 256
                fillerElements = 256 - mod(headerLength + 4, 256);
                for i = 1:fillerElements
                    fprintf(fid, '\n');
                end;
            end;
            
            fwrite(fid, '##}');
            fprintf(fid, '\n');
            
            fwrite(fid, convertedStack(:), 'uint16', machineFormat);
            
            fclose(fid);
        else
            writeImage(convertedStack, outputName);
        end;
    end;
end;

if matlabpool('size') > 0
    matlabpool('close');
end;