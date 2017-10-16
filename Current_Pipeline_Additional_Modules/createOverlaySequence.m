background = 100;

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..') || strcmp(candidates(i).name, 'Overlays')
        candidates(i) = [];
    else    
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb'], 'file') ~= 2;
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
    
    frame = single(readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_xy.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_xy.klb']));
    frame = frame ./ 4;
    frames{i, 1} = uint16(frame);
    sizeArray(i, 1:2) = [size(frame, 1) size(frame, 2)];
    
    frame = single(readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.affine.trsf_yz.klb']));
    frame = frame + single(readImage([candidates(i).name '\SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.affine.trsf_yz.klb']));
    frame = frame ./ 4;
    frames{i, 2} = uint16(frame);
    sizeArray(i, 3:4) = [size(frame, 1) size(frame, 2)];
end;

xMax = max(sizeArray(:, 1));
yMax = max(sizeArray(:, 2));

xyOverlays = ones(xMax, yMax, numel(candidates), 'uint16') .* background;

for i = 1:numel(candidates)
    xyOverlays(...
        (round((xMax - sizeArray(i, 1)) / 2) + 1):(round((xMax - sizeArray(i, 1)) / 2) + sizeArray(i, 1)), ...
        (round((yMax - sizeArray(i, 2)) / 2) + 1):(round((yMax - sizeArray(i, 2)) / 2) + sizeArray(i, 2)), ...
        i) = frames{i, 1};
end;

writeImage(xyOverlays, 'Overlays\xyOverlays.klb');

yMax = max(sizeArray(:, 3));
zMax = max(sizeArray(:, 4));

yzOverlays = ones(yMax, zMax, numel(candidates), 'uint16') .* background;

for i = 1:numel(candidates)
    yzOverlays(...
        (round((yMax - sizeArray(i, 3)) / 2) + 1):(round((yMax - sizeArray(i, 3)) / 2) + sizeArray(i, 3)), ...
        (round((zMax - sizeArray(i, 4)) / 2) + 1):(round((zMax - sizeArray(i, 4)) / 2) + sizeArray(i, 4)), ...
        i) = frames{i, 2};
end;

writeImage(yzOverlays, 'Overlays\yzOverlays.klb');