%% general parameters

% shiftsName    = 'Shift Data' filesep 'shifts_160113.mat';
% inputString2  = '_timeFused_wavelet' filesep 'SPM00_TM';
% inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-01-13' filesep 'Mmu_E1_GaleGFP_20160113_171610.corrected' filesep 'Results' filesep 'TimeFused.Wavelet' filesep 'Mmu_E1_GaleGFP.TM';
% allTimepoints = 0:322;
% timeSelection = 0:322;
% channels      = 0:1;
% copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-01-13 Temporal Registration Assessment' filesep 'TimeFused.Wavelet.Shifted.Projections';

% shiftsName    = 'Shift Data' filesep 'shifts_160614.mat';
% inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-06-14' filesep 'Mmu_E1_Foxa2eGFPxTmCherry_01_20160614_121243.corrected' filesep 'Results' filesep 'TimeFused.Geometric' filesep 'Mmu_E1_Foxa2eGFPxTmCherry.TM';
% inputString2  = '_timeFused_blending' filesep 'SPM00_TM';
% allTimepoints = 0:258;
% timeSelection = 0:258;
% channels      = 0:1;
% copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-06-14 Temporal Registration Assessment' filesep 'TimeFused.Geometric.Shifted.Projections';

% shiftsName    = 'Shift Data' filesep 'shifts_160716.mat';
% inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-07-16' filesep 'Mmu_E1_GaleGFPxH2BmCherry_20160716_131448.corrected' filesep 'Results' filesep 'TimeFused.Blending' filesep 'Mmu_E1_GaleGFPxH2BmCherry.TM';
% inputString2  = '_timeFused_blending' filesep 'SPM00_TM';
% allTimepoints = 0:271;
% timeSelection = 0:271;
% channels      = 0:1;
% copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-07-16 Temporal Registration Assessment' filesep 'TimeFused.Blending.Shifted.Projections';

% shiftsName    = 'Shift Data' filesep 'shifts_160716.mat';
% inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-07-16' filesep 'Mmu_E1_GaleGFPxH2BmCherry_20160716_131448.corrected' filesep 'Results' filesep 'TimeFused.Wavelet' filesep 'Mmu_E1_GaleGFPxH2BmCherry.TM';
% inputString2  = '_timeFused_blending' filesep 'SPM00_TM';
% allTimepoints = 0:271;
% timeSelection = 0:271;
% channels      = 0:1;
% copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-07-16 Temporal Registration Assessment' filesep 'TimeFused.Wavelet.Shifted.Projections';

% shiftsName    = 'Shift Data' filesep 'shifts_161215.mat';
% inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-12-15' filesep 'Mmu_E1_Fucci2_01_20161215_185256.corrected' filesep 'Results' filesep 'TimeFused.Geometric' filesep 'Mmu_E1_Fucci2.TM';
% inputString2  = '_timeFused_blending' filesep 'SPM00_TM';
% allTimepoints = 0:192;
% timeSelection = 0:192;
% channels      = 0:1;
% copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-12-15 Temporal Registration Assessment' filesep 'TimeFused.Geometric.Shifted.Projections';

shiftsName    = 'Shift Data' filesep 'shifts_161215.mat';
inputString1  = 'R:' filesep 'SV1' filesep 'KM_16-12-15' filesep 'Mmu_E1_Fucci2_01_20161215_185256.corrected' filesep 'Results' filesep 'TimeFused.Wavelet' filesep 'Mmu_E1_Fucci2.TM';
inputString2  = '_timeFused_wavelet' filesep 'SPM00_TM';
allTimepoints = 0:192;
timeSelection = 0:192;
channels      = 0:1;
copyLocation  = 'E:' filesep 'Mouse Development' filesep 'Time-Lapse Registration' filesep '16-12-15 Temporal Registration Assessment' filesep 'TimeFused.Wavelet.Shifted.Projections';

inputString3  = '_CM00_CM01_CHN';
inputString4  = '.fusedStack.klb';              % assumption in code below: stack name ends on 'corrected.klb'
projString1   = '.fusedStack_xyProjection.klb'; % assumption in code below: projection name ends on 'corrected.klb'
projString2   = '.fusedStack_xzProjection.klb'; % assumption in code below: projection name ends on 'corrected.klb'
projString3   = '.fusedStack_yzProjection.klb'; % assumption in code below: projection name ends on 'corrected.klb'

sourceType    = 2; % 0 for manual shift data (pairwise relative vectors)
                   % 1 for automated shift data (global vectors) computed by Gregoire's code
                   % 2 for automated shift data (global vectors) computed by Leo's cross-correlation
makeLocalCopy = 1; % put copy of projections in local folder

parallelJobs  = 12;

%% main loop

if makeLocalCopy && exist(copyLocation, 'dir') ~= 7
    mkdir(copyLocation);
end;

load(shiftsName);

if sourceType == 0
    shiftsAll = zeros(numel(allTimepoints), 4);
    shiftsAll(:, 1) = allTimepoints';
    for t = 1:numel(allTimepoints)
        shiftsAll(t, 2) = sum(shifts(shifts(:, 1) <= allTimepoints(t), 2));
        shiftsAll(t, 3) = sum(shifts(shifts(:, 1) <= allTimepoints(t), 3));
        shiftsAll(t, 4) = sum(shifts(shifts(:, 1) <= allTimepoints(t), 4));
    end;
    save(shiftsName, 'shifts', 'shiftsAll');
elseif sourceType == 1
    shiftsAll = zeros(numel(allTimepoints), 4);
    shiftsAll(:, 1) = allTimepoints';
    for t = 1:numel(allTimepoints)
        shiftsAll(t, 2) = round(shifts(shifts(:, 1) == allTimepoints(t), 2));
        shiftsAll(t, 3) = round(shifts(shifts(:, 1) == allTimepoints(t), 3));
        shiftsAll(t, 4) = round(shifts(shifts(:, 1) == allTimepoints(t), 4));
    end;
    save(shiftsName, 'shifts', 'shiftsAll');
else
    shiftsAll = zeros(numel(allTimepoints), 4);
    shiftsAll(:, 1) = allTimepoints';
    for t = 1:numel(allTimepoints)
        shiftsAll(t, 2) = -round(shifts(shifts(:, 1) == allTimepoints(t), 2));
        shiftsAll(t, 3) = -round(shifts(shifts(:, 1) == allTimepoints(t), 3));
        shiftsAll(t, 4) = -round(shifts(shifts(:, 1) == allTimepoints(t), 4));
    end;
    save(shiftsName, 'shifts', 'shiftsAll');
end;

minShiftX = min(shiftsAll(:, 2));
maxShiftX = max(shiftsAll(:, 2));
minShiftY = min(shiftsAll(:, 3));
maxShiftY = max(shiftsAll(:, 3));
minShiftZ = min(shiftsAll(:, 4));
maxShiftZ = max(shiftsAll(:, 4));

if matlabpool('size') > 0
    matlabpool('close');
end;
matlabpool(parallelJobs);

parfor t = 1:numel(timeSelection)
    for c = 1:numel(channels)
        disp(['Shifting image stack at TM' num2str(timeSelection(t)) ', channel ' num2str(channels(c))]);
        
        currentShift = shiftsAll(find(shiftsAll(:, 1) == timeSelection(t), 1), 2:4);
        timeString = num2str(timeSelection(t), '%.6d');
        channelString = num2str(channels(c), '%.2d');
        stackName = [inputString1 timeString inputString2 timeString inputString3 channelString inputString4];
        stack = readImage(stackName);
        
        shiftedStack = zeros(...
            size(stack, 1)+maxShiftX-minShiftX, ...
            size(stack, 2)+maxShiftY-minShiftY, ...
            size(stack, 3)+maxShiftZ-minShiftZ, 'uint16');
        
        shiftedStack(...
            (maxShiftX-currentShift(1)+1):(maxShiftX-currentShift(1)+size(stack, 1)), ...
            (maxShiftY-currentShift(2)+1):(maxShiftY-currentShift(2)+size(stack, 2)), ...
            (maxShiftZ-currentShift(3)+1):(maxShiftZ-currentShift(3)+size(stack, 3))) = stack;
        
        shiftedStackXYProj = max(shiftedStack, [], 3);
        shiftedStackXZProj = squeeze(max(shiftedStack, [], 2));
        shiftedStackYZProj = squeeze(max(shiftedStack, [], 1));
        
        shiftedStackName = [stackName(1:(end-3)) 'shifted.klb'];
        shiftedStackXYProjName = [stackName(1:(end-4)) '_xyProjection.shifted.klb'];
        shiftedStackXZProjName = [stackName(1:(end-4)) '_xzProjection.shifted.klb'];
        shiftedStackYZProjName = [stackName(1:(end-4)) '_yzProjection.shifted.klb'];
        writeImage(shiftedStack, shiftedStackName);
        writeImage(shiftedStackXYProj, shiftedStackXYProjName);
        writeImage(shiftedStackXZProj, shiftedStackXZProjName);
        writeImage(shiftedStackYZProj, shiftedStackYZProjName);
        
        if makeLocalCopy
            writeImage(shiftedStackXYProj, [copyLocation inputString2(20:end) timeString inputString3 channelString inputString4(1:(end-4)) '_xyProjection.shifted.klb']);
            writeImage(shiftedStackXZProj, [copyLocation inputString2(20:end) timeString inputString3 channelString inputString4(1:(end-4)) '_xzProjection.shifted.klb']);
            writeImage(shiftedStackYZProj, [copyLocation inputString2(20:end) timeString inputString3 channelString inputString4(1:(end-4)) '_yzProjection.shifted.klb']);
        end;
    end;
end;

if matlabpool('size') > 0
    matlabpool('close');
end;