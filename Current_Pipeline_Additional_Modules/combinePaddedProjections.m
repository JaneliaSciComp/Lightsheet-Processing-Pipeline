%% configuration

dataFolder = ['R:' filesep 'SV1' filesep 'KM_16-07-26' filesep 'Mmu_E1_GaleGFPxmKate2nls_01_20160726_162110.corrected' filesep 'Results' filesep 'TimeFused.Corrected'];
dirHeader  = 'Mmu_E1_GaleGFPxmKate2nls.TM';
dirFooter  = '_timeFused_blending';
fileHeader = 'SPM00_TM';
fileConfig = '_CM00_CM01_CHN01.fusedStack';
fileFooter = '.corrected.klb';

timepoints = 0:257;

outputDir  = ['E:' filesep 'Mouse Development' filesep 'Temp' filesep 'CHN01'];
dataTypes  = [1 0 1];
background = [0 100];
fullSearch = 1;

% dataFolder = 'R:' filesep 'SV1' filesep 'KM_16-01-13' filesep 'Mmu_E1_GaleGFP_20160113_171610.corrected' filesep 'Results' filesep 'TimeFused.Wavelet';
% dirHeader  = 'Mmu_E1_GaleGFP.TM';
% dirFooter  = '_timeFused_wavelet';
% fileHeader = 'SPM00_TM';
% fileConfig = '_CM00_CM01_CHN00.fusedStack';
% fileFooter = '.klb';
% 
% timepoints = 0:322;
% 
% outputDir  = 'E:' filesep 'Mouse Development' filesep 'Temp';
% dataTypes  = [1 1 1];
% background = [0 100];
% fullSearch = 1;

%% main loop

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end;

xyHeader = readKLBheader(...
    [dataFolder filesep '' dirHeader num2str(timepoints(end), '%.6d') dirFooter filesep '' ...
    fileHeader num2str(timepoints(end), '%.6d') fileConfig '_xyProjection' fileFooter]);
xzHeader = readKLBheader(...
    [dataFolder filesep '' dirHeader num2str(timepoints(end), '%.6d') dirFooter filesep '' ...
    fileHeader num2str(timepoints(end), '%.6d') fileConfig '_xzProjection' fileFooter]);
xSize = xyHeader.xyzct(1);
ySize = xyHeader.xyzct(2);
zSize = xzHeader.xyzct(2);

if fullSearch
    for t = 1:numel(timepoints)
        xzHeader = readKLBheader(...
            [dataFolder filesep '' dirHeader num2str(timepoints(t), '%.6d') dirFooter filesep '' ...
            fileHeader num2str(timepoints(t), '%.6d') fileConfig '_xzProjection' fileFooter]);
        zSize = max(zSize, xzHeader.xyzct(2));
    end;
end;

if background(1)
    xySequence = ones(xSize, ySize, numel(timepoints), 'uint16').*background(2);
    xzSequence = ones(xSize, zSize, numel(timepoints), 'uint16').*background(2);
    yzSequence = ones(ySize, zSize, numel(timepoints), 'uint16').*background(2);
else
    xySequence = zeros(xSize, ySize, numel(timepoints), 'uint16');
    xzSequence = zeros(xSize, zSize, numel(timepoints), 'uint16');
    yzSequence = zeros(ySize, zSize, numel(timepoints), 'uint16');
end;

disp(' ');

for t = 1:numel(timepoints)
    disp(['* Integrating time point ' num2str(timepoints(t))]);
    
    xyProjectionFilename = [dataFolder filesep '' dirHeader num2str(timepoints(t), '%.6d') dirFooter filesep '' ...
        fileHeader num2str(timepoints(t), '%.6d') fileConfig '_xyProjection' fileFooter];
    xzProjectionFilename = [dataFolder filesep '' dirHeader num2str(timepoints(t), '%.6d') dirFooter filesep '' ...
        fileHeader num2str(timepoints(t), '%.6d') fileConfig '_xzProjection' fileFooter];
    yzProjectionFilename = [dataFolder filesep '' dirHeader num2str(timepoints(t), '%.6d') dirFooter filesep '' ...
        fileHeader num2str(timepoints(t), '%.6d') fileConfig '_yzProjection' fileFooter];
    if dataTypes(1)
        xyProjection = readImage(xyProjectionFilename);
        xySequence(:, :, t) = xyProjection;
    end;
    if dataTypes(2)
        xzProjection = readImage(xzProjectionFilename);
        xzSequence(:, (1+floor((zSize-size(xzProjection, 2))/2)):(1+floor((zSize-size(xzProjection, 2))/2)+size(xzProjection, 2)-1), t) = xzProjection;
    end;
    if dataTypes(3)
        yzProjection = readImage(yzProjectionFilename);
        yzSequence(:, (1+floor((zSize-size(yzProjection, 2))/2)):(1+floor((zSize-size(yzProjection, 2))/2)+size(yzProjection, 2)-1), t) = yzProjection;
    end;
end;

disp(' ');
disp('Saving assembled image sequences');

if dataTypes(1)
    xySequenceFilename = [outputDir filesep 'xyProjections.klb'];
    writeImage(xySequence, xySequenceFilename);
end;
if dataTypes(2)
    xzSequenceFilename = [outputDir filesep 'xzProjections.klb'];
    writeImage(xzSequence, xzSequenceFilename);
end;
if dataTypes(3)
    yzSequenceFilename = [outputDir filesep 'yzProjections.klb'];
    writeImage(yzSequence, yzSequenceFilename);
end;

disp(' ');