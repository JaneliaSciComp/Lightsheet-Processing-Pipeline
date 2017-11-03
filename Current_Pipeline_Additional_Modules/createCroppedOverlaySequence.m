candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..') || strcmp(candidates(i).name, 'Overlays')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.xy.cropped.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.yz.cropped.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        end;
    end;
end;

if exist('Overlays', 'dir') ~= 7
    mkdir('Overlays');
end;

frames = cell(numel(candidates), 2);
sizeArray = zeros(numel(candidates), 4);

for i = 1:numel(candidates)
    timepoint = str2num(candidates(i).name(3:end));
    
    frames{i, 1} = readImage([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.xy.cropped.klb']);
    sizeArray(i, 1:2) = [size(frames{i, 1}, 1) size(frames{i, 1}, 2)];
    
    frames{i, 2} = readImage([candidates(i).name filesep 'TM' num2str(timepoint, '%.6d') '_overlay.yz.cropped.klb']);
    sizeArray(i, 3:4) = [size(frames{i, 2}, 1) size(frames{i, 2}, 2)];
end;

xMax = max(sizeArray(:, 1));
yMax = max(sizeArray(:, 2));

xyOverlays = zeros(xMax, yMax, numel(candidates), 'uint16');

for i = 1:numel(candidates)
    xyOverlays(...
        (round((xMax - sizeArray(i, 1)) / 2) + 1):(round((xMax - sizeArray(i, 1)) / 2) + sizeArray(i, 1)), ...
        (round((yMax - sizeArray(i, 2)) / 2) + 1):(round((yMax - sizeArray(i, 2)) / 2) + sizeArray(i, 2)), ...
        i) = frames{i, 1};
end;

writeImage(xyOverlays, 'Overlays' filesep 'xyCroppedOverlays.klb');

yMax = max(sizeArray(:, 3));
zMax = max(sizeArray(:, 4));

yzOverlays = zeros(yMax, zMax, numel(candidates), 'uint16');

for i = 1:numel(candidates)
    yzOverlays(...
        (round((yMax - sizeArray(i, 3)) / 2) + 1):(round((yMax - sizeArray(i, 3)) / 2) + sizeArray(i, 3)), ...
        (round((zMax - sizeArray(i, 4)) / 2) + 1):(round((zMax - sizeArray(i, 4)) / 2) + sizeArray(i, 4)), ...
        i) = frames{i, 2};
end;

writeImage(yzOverlays, 'Overlays' filesep 'yzCroppedOverlays.klb');